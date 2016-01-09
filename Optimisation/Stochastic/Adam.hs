{-# LANGUAGE DeriveGeneric #-}
{-# LANGUAGE MultiParamTypeClasses #-}

module Optimisation.Stochastic.Adam (
    adam,         -- the optimisation function
    AdamParams(..)
  ) where

import Control.Monad
import Control.Monad.Trans.Class
import Control.Monad.Trans.State.Strict
import Control.Monad.Trans.Reader
import Data.Serialize (Serialize)
import Data.Default
import GHC.Generics
import System.Random

import Optimisation.Stochastic.Common

-- Adam algorithm for stochastic functions optimisation

-- Per dimension needed info
data Dim = Dim {
               crt      :: !Double,     -- current value
               del      :: !Double,     -- gradient calculation delta
               gra      :: !Double,     -- current gradient value
               nxt      :: !Double,     -- next value
               frm      :: !Double,     -- first moment estimate
               srm      :: !Double      -- second raw moment estimate
           }
           deriving (Show, Generic)

instance Dimension Dim where
    dimCrt = crt
    dimDlt = del
    setGrd g dim = dim { gra = g }

data GState = GState {
                 save :: Maybe FilePath, -- file to use for checkpoints
                 rfil :: Maybe FilePath, -- running file name (delete to stop asap)
                 dims :: [Dim],          -- per dimension state
                 step :: !Int,           -- current step
                 sstp :: !Int,           -- steps close to termination
                 cnxx :: !Double,        -- distance between crt & nxt
                 b1n, b2n :: !Double     -- powers of decays
             } deriving Generic

instance Default GState where
    def = GState { save = Nothing, dims = [], step = 1, sstp = 0, cnxx = 0, b1n = 1, b2n = 1,
                   rfil = Nothing }

-- Some parameter of the algorithm (remain constant during one execution):
data AdamParams = AdamParams {
                  verb :: !Bool,       -- verbose?
                  cent :: !Bool,       -- central second momentum?
                  alfa :: !Double,     -- alpha step size
                  bet1 :: !Double,     -- exponential decay for first order momentum
                  bet2 :: !Double,     -- exponential decay for second order momentum
                  epsi :: !Double,     -- to limit the update step
                  xstp :: !Double,     -- termination when max(n) (xn - xc) < xstp
                  ngrs :: !Int,        -- number of gradients per step
                  nmax :: !Int,        -- maximum number of steps at all
                  nstp :: !Int         -- steps to stay close for termination
              }
              deriving (Show, Generic)

-- These instances are used to save/load state
instance Serialize Dim
instance Serialize AdamParams
instance Serialize GState

-- Default params
instance Default AdamParams where
    def = AdamParams { verb = False, alfa = 0.1, bet1 = 0.9, bet2 = 0.999, epsi = 1,
                       xstp = 1/1000, nmax = 1000, nstp = 3, cent = False, ngrs = 1 }

instance OptimState GState AdamParams where
    iniCurr dlucs s = s { dims = map (uncurry startDim) dlucs }
    getCurr       s = map crt (dims s)
    setSave file  s = s { save = Just file }
    getSave       s = save s
    setRun  file  s = s { rfil = Just file }
    getRun        s = rfil s

-- Initialisation per dimension
-- Expects: current (start point), gradient delta
startDim :: Double -> Double -> Dim
startDim c d = Dim { crt = c, del = d, gra = 0, nxt = 0, frm = 0, srm = 0 }

nextMoment :: AdamParams -> Double -> Dim -> Dim
nextMoment pars b1n dim = dim { frm = m1, srm = m2 }
    where m1 = bet1 pars * frm dim + (1 - bet1 pars) * gra dim
          m2 | cent pars = let m' = m1 / (1 - b1n)
                               gd = gra dim - m'
                           in bet2 pars * srm dim + (1 - bet2 pars) * gd * gd
             | otherwise = bet2 pars * srm dim + (1 - bet2 pars) * gra dim * gra dim

-- Adding here means: we maximise
nextPoint :: AdamParams -> Double -> Double -> Dim -> Dim
nextPoint pars b1n b2n dim = dim { nxt = x }
    where m = frm dim / (1 - b1n)
          v = srm dim / (1 - b2n)
          x = crt dim + alfa pars * m / (sqrt v + epsi pars)   -- v = 0 ==> step = alpha * m

-- Prepare next step, updating and projecting crt point
nextStep :: Dim -> Dim
nextStep dim = dim { crt = nxt dim }

-- The norm for the termination condition
dist1 :: [Double] -> [Double] -> Double
dist1 as bs = maximum $ map abs $ zipWith subtract as bs

adam :: Stochastic -> Bool -> Maybe AdamParams -> StartOptions -> IO [Double]
adam play rf mparams staopts = optimise play optimStep rf mparams staopts

optimStep :: Stochastic -> Optim GState AdamParams Bool
optimStep play = checkRunning $ do
    stat <- lift get
    params <- ask
    checkEnd (verb params) (step stat > nmax params) "Termination: nmax steps reached" $
      checkEnd (verb params) (cnxx stat < xstp params && sstp stat >= nstp params)
        "Termination: optimum near enough" $ do
           --  info $ "Crt = " ++ show (map crt (dims stat))
           --  info $ "Dims = " ++ show (dims stat)
           when (verb params) $ do
               info ""
               info $ "*** Step " ++ show (step stat) ++ " ***"
               info $ "Crt = " ++ show (map crt $ dims stat)
               info $ "Calculate gradient..."
           -- Calculate gradient at current point
           dgs <- lift . lift $ calcNGrad (verb params) play (ngrs params) (dims stat)
           when (verb params) $ info $ "Gra = " ++ show (map gra dgs)
           let -- update the momentum per dimension
               b1n' = b1n stat * bet1 params
               dms = map (nextMoment params b1n') dgs
           when (verb params) $ do
               info $ "Mo1 = " ++ show (map frm dms)
               info $ "Mo2 = " ++ show (map srm dms)
           let -- next point
               b2n' = b2n stat * bet2 params
               ds = map (nextPoint params b1n' b2n') dms
           when (verb params) $ info $ "Nxt = " ++ show (map nxt ds)
           let -- next point, distance in parameter space (for termination)
               xx  = dist1 (map crt ds) (map nxt ds)
               dns = map nextStep ds
               near | xx < xstp params = sstp stat + 1
                    | otherwise        = 0
           when (verb params) $ info $ "Dist1 = " ++ show xx
           lift $ put stat { dims = dns, step = step stat + 1, sstp = near,
                             cnxx = xx, b1n = b1n', b2n = b2n' }
           return False

------------------------
-- Test optimisations --
------------------------

-- 1 dimension, no noise: a simple quadratic funtion with maximum in -10
oneDim :: [Double] -> IO Double
oneDim (x:_) = return $ 1 - 0.1 * (x+10) * (x+10)

maxOneDimExact n = adam oneDim False
                        (Just def { verb = True, alfa = 1, nmax = n, xstp = 1E-6 })
                        (SOStartNoCheckpoint [(50, 1)])

-- Same with noise: noise with a general level
withNoise :: Double -> Stochastic -> [Double] -> IO Double
withNoise noise play xs = do
    y <- play xs
    r <- getStdRandom (randomR (-noise, noise))
    return $ y + r

maxOneDimNoise n = adam (withNoise 0.5 oneDim) False
                        (Just def { verb = True, alfa = 1, nmax = n, xstp = 1E-6 })
                        (SOStartNoCheckpoint [(0, 1)])

-- 2 dimension: negated banana, no noise, (global) maximum at (1, 1)
banana :: [Double] -> IO Double
banana (x:y:_) = return $ negate $ x1 * x1 + 10 * y1 * y1
    where x1 = 1 - x
          x2 = x * x
          y1 = y - x2

maxBanana n k dx v = adam banana False
                         (Just def { verb = v, ngrs = k, nmax = n, xstp = 1E-6 })
                         (SOStartNoCheckpoint [(-5, dx), (-5, dx)])
