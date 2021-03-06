{-# LANGUAGE DeriveGeneric #-}

module Radexp (
    radexp,         -- the optimisation function
    Params(..),
    StartOptions(..),
    defRadexpParams,
  ) where

import Control.Monad
import Control.Monad.Trans.Class
import Control.Monad.Trans.State.Strict
import qualified Data.ByteString.Lazy as BL
import Data.Serialize (Serialize, encodeLazy, decodeLazy)
import GHC.Generics
import System.Directory
import System.FilePath
import System.Random

-- Radexp algorithm for stochastic functions optimisation -- own invention!

-- We need to let some external process to get a sample of the function
-- at some position, and this will happen generally in IO
type Stochastic = [Double] -> IO Double
type ORange = (Double, Double)
data StartOptions
    = SOStartWithCheckpoint FilePath [ORange]
    | SOStartNoCheckpoint            [ORange]
    | SORestart             FilePath

-- Per dimension needed info
data Dim = Dim {
               prv      :: !Double,     -- previous value
               crt      :: !Double,     -- current value
               del      :: !Double,     -- gradient calculation delta
               grm      :: !Double,     -- gradient momentum
               dex      :: !Double,     -- delta x
               nxt      :: !Double,     -- next value
               hes      :: !Double      -- hessian element
           }
           deriving (Show, Generic)

data GState = GState {
                 pars :: Params,         -- we include the params in the state record
                 save :: Maybe FilePath, -- file to use for checkpoints
                 rfil :: FilePath,       -- running file name (delete to stop asap)
                 dims :: [Dim],          -- per dimension state
                 step :: !Int,           -- current step
                 sstp :: !Int,           -- steps close to termination
                 cnxx :: !Double         -- distance between crt & nxt
             } deriving Generic

-- Our monad stack
type Optim a = StateT GState IO a

-- Some parameter of the algorithm (remain constant during one execution):
data Params = Params {
                  verb :: !Bool,       -- verbose?
                  beta :: !Double,     -- exponential decay
                  gama :: !Double,     -- limit target per dimension
                  epsi :: !Double,     -- to limit the update step
                  xstp :: !Double,     -- termination when max(n) (xn - xc) < xstp
                  nmax :: !Int,        -- maximum number of steps at all
                  nstp :: !Int         -- steps to stay close for termination
              }
              deriving (Show, Generic)

-- These instances are used to save/load state
instance Serialize Dim
instance Serialize Params
instance Serialize GState

-- Default params
defRadexpParams :: Params
defRadexpParams = Params { verb = False, beta = 0.9, gama = 1, epsi = 1E-8, xstp = 1E-8,
                           nmax = 1000, nstp = 3 }

-- Initialisation per dimension
-- Expects: current (start point), gradient delta
startDim :: Double -> Double -> Dim
startDim c d = Dim { prv = c, crt = c, del = d, grm = 0, dex = 0, nxt = 0, hes = 0 }

-- r has same signum as gradiend, so adding here means, we maximise
nextPoint :: Params -> Dim -> Dim
nextPoint pars dim = dim { dex = dx, nxt = x }
    where r  = grm dim / (abs (hes dim) + epsi pars)
          dx = gama pars * r / (1 + abs r)    -- we don't trust big steps
          x  = crt dim + dx

-- Prepare next step, updating and projecting crt point
nextStep :: Dim -> Dim
nextStep dim = dim { prv = crt dim, crt = nxt dim }

-- The norm for the termination condition
dist1 :: [Double] -> [Double] -> Double
dist1 as bs = maximum $ map abs $ zipWith subtract as bs

-- Uniformly random +/-1
randSign :: IO Int
randSign = do
    r <- getStdRandom (randomR (1,2))
    case r of
        1 -> return r
        2 -> return (-1)

-- Calculate gradient and "hessian"
calcGrad :: Params -> Stochastic -> [Dim] -> IO [Dim]
calcGrad pars play ds = do
    let xs = map crt ds
    dx <- forM ds $ \dim -> do
              d <- randSign
              return $ fromIntegral d * del dim
    let xp = zipWith (+)      dx xs
        xm = zipWith subtract dx xs
    fp <- play xp
    when (verb pars) $ putStrLn $ "Func + at " ++ show xp ++ ": " ++ show fp
    fm <- play xm
    when (verb pars) $ putStrLn $ "Func - at " ++ show xm ++ ": " ++ show fm
    let cgr = map (grad (fp - fm)) $ zipWith subtract xm xp     -- current gradient
        pxs = map prv ds
        pgr = map grm ds     -- grm is considered as previous gradient, at the previous point
        mgr = zipWith (\p c -> p * beta pars + c * (1 - beta pars)) pgr cgr
        hes = zipWith grad (zipWith subtract pgr mgr) $ zipWith subtract pxs xs
        dgs = zipWith setgh ds $ zip mgr hes
    return dgs
    where grad df dx
              -- | abs dx < epsi pars = 0   -- is this ok? It ignores noise, which is good
              | abs dx == 0 = 0
              | otherwise   = df / dx
          setgh dim (g, h) = dim { grm = g, hes = h }

-- Low limit of denumerators (to avoid devide by 0 or too high results)
limEpsi :: Double -> Double -> Double
limEpsi eps x | x > eps || x < -eps = x
              | otherwise           = eps

doUntil :: Monad m => m Bool -> m ()
doUntil act = go
    where go = do
              r <- act
              if r then return () else go

radexp :: Stochastic -> Maybe Params -> StartOptions -> IO [Double]
radexp play mparams staopts = do
    stat <- case staopts of
                SOStartWithCheckpoint file dlucs -> do
                    let params = maybe defRadexpParams id mparams
                        ds = map (uncurry startDim) dlucs
                        stat = GState { pars = params, save = Just file, dims = ds,
                                        step = 1, sstp = 0, cnxx = 0, rfil = ""
                                      }
                    return stat
                SOStartNoCheckpoint        dlucs -> do
                    let params = maybe defRadexpParams id mparams
                        ds = map (uncurry startDim) dlucs
                        stat = GState { pars = params, save = Nothing, dims = ds,
                                        step = 1, sstp = 0, cnxx = 0, rfil = ""
                                      }
                    return stat
                SORestart             file       -> do
                    stat' <- restorePoint file
                    let params = maybe (pars stat') id mparams  -- we can override params
                        stat = stat' { pars = params, save = Just file }
                    return stat
    -- Write an empty file with a random number in name to be used for stopping
    -- an optimisation run as soon as possible, but in a smooth way
    -- If the file will be deleted, the optimisation will stop after the current step ends
    -- If you have no checkpoint file, you will not be able to continue from there!
    r <- runFileNum
    let rfile = "running-" ++ show r
    BL.writeFile rfile BL.empty
    (ds, fs) <- runStateT (spsa play) stat { rfil = rfile }
    return ds

runFileNum :: IO Int
runFileNum = getStdRandom (randomR (1000, 9999))

spsa :: Stochastic -> Optim [Double]
spsa play = do
    doUntil $ optimStep play
    gets $ map crt . dims

optimStep :: Stochastic -> Optim Bool
optimStep play = checkRunning $ do
    stat <- get
    let params = pars stat
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
           let nn = fromIntegral $ step stat
           -- Calculate gradient at current point
           dgs <- lift $ calcGrad params play (dims stat)
           when (verb params) $ do
               info $ "Prv = " ++ show (map prv dgs)
               info $ "Grm = " ++ show (map grm dgs)
               info $ "Hes = " ++ show (map hes dgs)
           let -- next point
               ds = map (nextPoint params) dgs
           when (verb params) $ do
               info $ "Dex = " ++ show (map dex ds)
               info $ "Nxt = " ++ show (map nxt ds)
           let -- next point, distance in parameter space (for termination)
               xx  = dist1 (map crt ds) (map nxt ds)
               dns = map nextStep ds
               near | xx < xstp params = sstp stat + 1
                    | otherwise        = 0
           when (verb params) $ info $ "Dist1 = " ++ show xx
           put stat { dims = dns, step = step stat + 1, sstp = near, cnxx = xx }
           checkPoint
           return False

info :: String -> Optim ()
info = lift . putStrLn

checkEnd :: Bool -> Bool -> String -> Optim Bool -> Optim Bool
checkEnd verb cond mes act
    | cond      = do
        when verb $ info mes
        return True
    | otherwise = act

checkPoint :: Optim ()
checkPoint = do
    s <- get
    case save s of
        Nothing -> return ()
        Just cp -> lift $ do
            let cpn = addExtension cp "new"
            BL.writeFile cpn $ encodeLazy s
            renameFile cpn cp

restorePoint :: FilePath -> IO GState
restorePoint file = do
    bs <- BL.readFile file
    case decodeLazy bs of
        Left estr -> error $ "Decode: " ++ estr
        Right sta -> return sta

checkRunning :: Optim Bool -> Optim Bool
checkRunning act = do
    s <- get
    rfe <- lift $ doesFileExist (rfil s)
    if rfe then act else do
       lift $ putStrLn $ "Running file " ++ rfil s ++ " not found, stop optimisation"
       -- put s { rfdl = False }
       return True

------------------------
-- Test optimisations --
------------------------

-- 1 dimension, no noise: a simple quadratic funtion with maximum in -10
oneDim :: [Double] -> IO Double
oneDim (x:_) = return $ 1 - 0.1 * (x+10) * (x+10)

maxOneDimExact n = radexp oneDim
                        (Just defRadexpParams { verb = True, nmax = n, xstp = 1E-6 })
                        (SOStartNoCheckpoint [(50, 1)])

-- Same with noise: noise with a general level
withNoise :: Double -> Stochastic -> [Double] -> IO Double
withNoise noise play xs = do
    y <- play xs
    r <- getStdRandom (randomR (-noise, noise))
    return $ y + r

maxOneDimNoise n = radexp (withNoise 0.5 oneDim)
                        (Just defRadexpParams { verb = True, nmax = n, xstp = 1E-6 })
                        (SOStartNoCheckpoint [(0, 1)])

-- 2 dimension: negated banana, no noise, (global) maximum at (1, 1)
banana :: [Double] -> IO Double
banana (x:y:_) = return $ negate $ x1 * x1 + 10 * y1 * y1
    where x1 = 1 - x
          x2 = x * x
          y1 = y - x2

maxBanana n = radexp banana
                   (Just defRadexpParams { nmax = n, xstp = 1E-6 })
                   (SOStartNoCheckpoint [(5, 1), (5, 1)])
