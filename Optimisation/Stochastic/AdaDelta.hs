{-# LANGUAGE DeriveGeneric #-}
{-# LANGUAGE MultiParamTypeClasses #-}

module Optimisation.Stochastic.AdaDelta (
    adaDelta,         -- the optimisation function
    AdaDeltaParams(..)
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

-- AdaDelta algorithm for stochastic functions optimisation

-- Per dimension needed info
data Dim = Dim {
               crt      :: !Double,     -- current value
               del      :: !Double,     -- gradient calculation delta
               gra      :: !Double,     -- current gradient value
               grm      :: !Double,     -- gradient momentum
               nxt      :: !Double,     -- next value
               dex      :: !Double,     -- last delta x
               edx2     :: !Double,     -- RMS of delta x
               eg2      :: !Double      -- RMS of gradient
           }
           deriving (Show, Generic)

instance Dimension Dim where
    dimCrt = crt
    dimDlt = del
    setGrd g dim = dim { gra = g }

data GState = GState {
                 save :: Maybe FilePath, -- file to use for checkpoints
                 rfil :: Maybe FilePath,       -- running file name (delete to stop asap)
                 dims :: [Dim],          -- per dimension state
                 step :: !Int,           -- current step
                 sstp :: !Int,           -- steps close to termination
                 cnxx :: !Double,        -- distance between crt & nxt
                 betn :: !Double         -- bias correction for exponetial decay vars
             } deriving Generic

instance Default GState where
    def = GState { save = Nothing, rfil = Nothing, dims = [], step = 1, sstp = 5,
                   cnxx = 1E-6, betn = 1 }

-- Some parameter of the algorithm (remain constant during one execution):
data AdaDeltaParams = AdaDeltaParams {
                  verb :: !Bool,       -- verbose?
                  alfa :: !Double,     -- step size factor
                  beta :: !Double,     -- exponential decay for RMS
                  epsi :: !Double,     -- to limit the update step
                  xstp :: !Double,     -- termination when max(n) (xn - xc) < xstp
                  ngrd :: !Int,        -- number of gradients per point
                  nmax :: !Int,        -- maximum number of steps at all
                  nstp :: !Int         -- steps to stay close for termination
              }
              deriving (Show, Generic)

-- These instances are used to save/load state
instance Serialize Dim
instance Serialize AdaDeltaParams
instance Serialize GState

-- Default params
instance Default AdaDeltaParams where
    def = AdaDeltaParams { verb = False, alfa = 1, beta = 0.95, epsi = 1E-6,
                           xstp = 1E-6, ngrd = 1, nmax = 1000, nstp = 10 }

instance OptimState GState AdaDeltaParams where
    iniCurr dlucs s = s { dims = map (uncurry startDim) dlucs }
    getCurr       s = map crt (dims s)
    setSave file  s = s { save = Just file }
    getSave       s = save s
    setRun  file  s = s { rfil = Just file }
    getRun        s = rfil s

-- Initialisation per dimension
-- Expects: current (start point), gradient delta
startDim :: Double -> Double -> Dim
startDim c d = Dim { crt = c, del = d, gra = 0, grm = 0, nxt = 0, dex = 0, edx2 = 0, eg2 = 0 }

rms :: AdaDeltaParams -> Double -> Double
rms pars x = sqrt (x + epsi pars)

corm :: Double -> Double -> Double
corm b x = x / (1 - b)

accumGrad :: AdaDeltaParams -> Double -> Dim -> Dim
accumGrad pars bn dim = dim { grm = g, eg2 = x }
    where g  = beta pars * grm dim + (1 - beta pars) * gra dim
          g' = corm bn g
          x  = beta pars * eg2 dim + (1 - beta pars) * g' * g'

compUpdate :: AdaDeltaParams -> Double -> Dim -> Dim
compUpdate pars bn dim = dim { dex = x }
    where g  = corm bn (grm dim)
          g2 = corm bn (eg2 dim)
          x2 = corm bn (edx2 dim)
          x  = alfa pars * g * rms pars x2 / rms pars g2

accumUpdate :: AdaDeltaParams -> Dim -> Dim
accumUpdate pars dim = dim { edx2 = x }
    where x = beta pars * edx2 dim + (1 - beta pars) * dex dim * dex dim

-- Adding here means: we maximise
nextPoint :: AdaDeltaParams -> Dim -> Dim
nextPoint pars dim = dim { nxt = x }
    where x = crt dim + dex dim

-- Prepare next step, updating and projecting crt point
nextStep :: Dim -> Dim
nextStep dim = dim { crt = nxt dim }

-- The norm for the termination condition
dist1 :: [Double] -> [Double] -> Double
dist1 as bs = maximum $ map abs $ zipWith subtract as bs

doUntil :: Monad m => m Bool -> m ()
doUntil act = go
    where go = do
              r <- act
              if r then return () else go

adaDelta :: Stochastic -> Bool -> Maybe AdaDeltaParams -> StartOptions -> IO [Double]
adaDelta play rf mparams staopts = optimise play optimStep rf mparams staopts

optimStep :: Stochastic -> Optim GState AdaDeltaParams Bool
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
           let bn = betn stat * beta params
           -- Calculate gradient at current point
           dgs <- lift . lift $ calcNGrad (verb params) play (ngrd params) (dims stat)
           when (verb params) $ info $ "Gra = " ++ show (map gra dgs)
           let -- accumulate gradients per dimension
               dag = map (accumGrad params bn) dgs
           when (verb params) $ do
               info $ "Grm = " ++ show (map grm dag)
               info $ "Eg2 = " ++ show (map eg2 dag)
           let -- compute updates per dimension
               dcu = map (compUpdate params bn) dag
           when (verb params) $ info $ "Dx  = " ++ show (map dex dcu)
           let -- accumulate updates per dimension
               dau = map (accumUpdate params) dcu
           when (verb params) $ info $ "Edx = " ++ show (map edx2 dau)
           let -- next point
               ds = map (nextPoint params) dau
           when (verb params) $ info $ "Nxt = " ++ show (map nxt ds)
           let -- next point, distance in parameter space (for termination)
               xx  = dist1 (map crt ds) (map nxt ds)
               dns = map nextStep ds
               near | xx < xstp params = sstp stat + 1
                    | otherwise        = 0
           when (verb params) $ info $ "Dist1 = " ++ show xx
           lift $ put stat { dims = dns, step = step stat + 1, sstp = near, cnxx = xx, betn = bn }
           return False

------------------------
-- Test optimisations --
------------------------

-- 1 dimension, no noise: a simple quadratic funtion with maximum in -10
oneDim :: [Double] -> IO Double
oneDim (x:_) = return $ 1 - 0.1 * (x+10) * (x+10)

maxOneDimExact n = adaDelta oneDim False
                        (Just def { verb = True, nmax = n, xstp = 1E-6 })
                        (SOStartNoCheckpoint [(50, 1)])

-- Same with noise: noise with a general level
withNoise :: Double -> Stochastic -> [Double] -> IO Double
withNoise noise play xs = do
    y <- play xs
    r <- getStdRandom (randomR (-noise, noise))
    return $ y + r

maxOneDimNoise n = adaDelta (withNoise 0.5 oneDim) False
                        (Just def { verb = True, nmax = n, xstp = 1E-6 })
                        (SOStartNoCheckpoint [(0, 1)])

-- 2 dimension: negated banana, no noise, (global) maximum at (1, 1)
banana :: [Double] -> IO Double
banana (x:y:_) = return $ negate $ x1 * x1 + 10 * y1 * y1
    where x1 = 1 - x
          x2 = x * x
          y1 = y - x2

maxBanana n = adaDelta banana False
                   (Just def { verb = True, nmax = n, xstp = 1E-6 })
                   (SOStartNoCheckpoint [(5, 1), (5, 1)])
