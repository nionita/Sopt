module SSSPSA (
    ssSPSA,         -- the optimisation function
    Params(..),
    defSpsaParams,
  ) where

import Control.Monad
import Control.Monad.Trans.Class
import Control.Monad.Trans.State.Strict
import System.Random

-- Scaled and shifted SPSA algorithm for stochastic functions optimisation
-- The SS-KW algorithm could also be implemented, but there are small differences
-- in treating the limit regions: KW is doing the gradient only in the posivite side
-- Will have to abstract this away somehow...

-- We need to let some external process to get a sample of the function
-- at some position, and this will happen generally in IO
type Stochastic = [Double] -> IO Double

-- Per dimension needed info
data Dim = Dim {
               low, hig :: !Double,     -- range of dimension
               alf, bet :: !Double,     -- parameters of a(n) sequence
               vak      :: !Double,     -- shift upper bound
               gam      :: !Double,     -- parameter of c(n) sequence
               cmx      :: !Double,     -- maximum for c(n) sequence
               crt      :: !Double,     -- current value
               gra      :: !Double,     -- current gradient value
               nxt      :: !Double,     -- next value
               ash      :: !Int,        -- number of a shifts so far
               csc      :: !Int         -- number of c sclaes so far
           }
           deriving Show

data GState = GState {
                 pars :: Params,         -- we include the params in the state record
                 dims :: [Dim],          -- per dimension state
                 oscs :: [Bool],         -- already oscillated dimensions
                 step :: !Int,           -- current step
                 grde :: !Int,           -- number of gradient estimations
                 sstp :: !Int,           -- steps close to termination
                 cnxx :: !Double         -- distance between crt & nxt
             }

-- Our monad stack
type Optim a = StateT GState IO a

-- Some parameter of the algorithm (remain constant during one execution):
data Params = Params {
                  verb :: !Bool,       -- verbose?
                  phia :: !Double,     -- maximul alpha scaling factor
                  c0   :: !Double,     -- define the maximum gamma: fraction from initial interval
                  gam0 :: !Double,     -- gamma scaling factor
                  gami :: !Double,     -- fraction from initial interval for initial gamma
                  xstp :: !Double,     -- termination when norm1 (xn -xc) < xstp
                  nmax :: !Int,        -- maximum number of steps at all
                  h0   :: !Int,        -- max steps for oscillations (a scaling)
                  ka   :: !Int,        -- maximum number of a(n) shifts per dimension
                  kc   :: !Int,        -- maximum number of c(n) scale ups per dimension
                  gmax :: !Int,        -- max allowed gradients for occilation phase
                  mmax :: !Int,        -- max step to which shifts/scalings are allowed
                  nstp :: !Int         -- steps to stay close for termination
              }
              deriving Show

-- Default params
defSpsaParams :: Params
defSpsaParams = Params { verb = False, phia = 10, c0 = 0.2, gam0 = 2, gami = 20, xstp = 1/1000000,
                         nmax = 1000, h0 = 4, ka = 50, kc = 50, gmax = 20, mmax = 1000, nstp = 3 }

-- Initialisation per dimension
-- Expects: low, high and current (start point)
startDim :: Params -> Double -> Double -> Double -> Dim
startDim par l h c
    = Dim { low = l, hig = h, alf = 1, bet = 0, vak = 10,
            gam = (h - l) / gami par, cmx = c0 par * (h - l), crt = c,
            gra = 0, nxt = 0, ash = 0, csc = 0 }

-- Calculate a(n) sequence
ak :: Dim -> Double -> Double
ak dim n = alf dim / (n + bet dim)

-- Calculate c(n) sequence
ck :: Dim -> Double -> Double
ck dim n = gam dim / (sqrt $ sqrt n)

-- Calculate lower bound
lbk :: Dim -> Double -> Double
lbk dim n = low dim + ck dim n

-- Calculate upper bound
ubk :: Dim -> Double -> Double
ubk dim n = hig dim - ck dim n

-- Alpha scaling step (per dimension)
ascalStep :: Params -> Double -> Dim -> (Bool, Dim)
ascalStep par n dim = (s, dim { alf = a })
    where (s, a) | gra dim > 0 = ascal par (ubk dim n) n dim
                 | gra dim < 0 = ascal par (lbk dim n) n dim
                 | otherwise   = (False, alf dim)     -- no scaling when gradient is 0

-- Calculating the alpha' scaling factor for alpha
-- Does not work for gradient 0
{-# INLINE ascal #-}
ascal :: Params -> Double -> Double -> Dim -> (Bool, Double)
ascal par lim n dim
    | f > phia par = (False, alf dim * phia par)  -- limited scaling
    | f < 1        = (True,  alf dim)             -- no scaling needed
    | otherwise    = (True,  al1)                 -- exact scaling
    where al1 = (lim - crt dim) / (nxt dim - crt dim)
          f   = al1 / alf dim

-- Beta shifting step (per dimension)
-- We expect here that the nxt field is not yet adjusted to the limits
bshiftStep :: Params -> Double -> Dim -> Dim
bshiftStep par n dim
    | ash dim >= ka par                 = dim
    | crt dim <= l1 && nxt dim > u2     = beshi u2
    | crt dim >= u1 && nxt dim < l2     = beshi l2
    | otherwise                         = dim
    where l1 = lbk dim n
          u1 = ubk dim n
          l2 = lbk dim (n+1)
          u2 = ubk dim (n+1)
          bshift delta = ceiling (alf dim * gra dim / delta - n - bet dim)
          beshi r = dim { bet = b, vak = v, ash = ash dim + 1 }
              where b' = fromIntegral $ bshift $ r - crt dim
                    (bp, v) | b' > vak dim = (vak dim, vak dim * 2)
                            | otherwise    = (b',      vak dim)
                    b = bet dim + bp

cscalStep :: Params -> Double -> Dim -> Dim
cscalStep par n dim
    | csc dim < kc par &&
      (  crt dim >= ubk dim n && nxt dim > crt dim
      || crt dim <= lbk dim n && nxt dim < crt dim) = dim { gam = gam dim * gap, csc = csc dim + 1 }
    | otherwise                                     = dim
    where gap = min (gam0 par) $ cmx dim / ck dim n

-- Adding here means: we maximise
nextPoint :: Double -> Dim -> Dim
nextPoint n dim = dim { nxt = x }
    where x = crt dim + ak dim n * gra dim

-- Prepare next step, updating and projecting crt point
nextStep :: Double -> Dim -> Dim
nextStep n dim = dim { crt = x }
    where l = lbk dim (n+1)
          u = ubk dim (n+1)
          x = min u $ max l $ nxt dim

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

-- Calculate gradient
calcGrad :: Bool -> Stochastic -> Double -> [Dim] -> IO [Dim]
calcGrad verb play n ds = do
    let xs = map crt ds
    dx <- forM ds $ \dim -> do
              d <- randSign
              return $ fromIntegral d * ck dim n
    let xp = zipWith (+)      dx xs
        xm = zipWith subtract dx xs
    fp <- play xp
    when verb $ putStrLn $ "Func + at " ++ show xp ++ ": " ++ show fp
    fm <- play xm
    when verb $ putStrLn $ "Func - at " ++ show xm ++ ": " ++ show fm
    let dgrad dim delta = dim { gra = (fp - fm) / delta }
        dgs = zipWith dgrad ds dx
    return dgs

doUntil :: Monad m => m Bool -> m ()
doUntil act = go
    where go = do
              r <- act
              if r then return () else go

ssSPSA :: Stochastic -> Params -> [((Double, Double), Double)] -> IO [Double]
ssSPSA play params dlucs = liftM fst $ runStateT (spsa play) stat
    where stat = GState { pars = params, dims = ds, oscs = os, step = 1,
                          grde = 0, sstp = 0, cnxx = 0 }
          ds = map (\((l, u), c) -> startDim params l u c) dlucs
          os = take (length ds) $ repeat False

spsa :: Stochastic -> Optim [Double]
spsa play = do
    doUntil $ scaleStep play
    doUntil $ shiftStep play
    gets $ map crt . dims

scaleStep :: Stochastic -> Optim Bool
scaleStep play = do
    stat <- get
    let params = pars stat
    if step stat > h0 params || grde stat > gmax params || all id (oscs stat)
       then do
           when (verb params) $ do
               if step stat > h0 params
                  then info $ "Scaling termination: h0 steps reached"
                  else if grde stat > gmax params
                          then info $ "Scaling termination: max gradient reached"
                          else info $ "Scaling termination: all dimensions scaled"
           return True
       else do
           when (verb params) $ do
               info ""
               info $ "Step " ++ show (step stat) ++ " (scale)"
               info $ "Crt = " ++ show (map crt $ dims stat)
               info $ "calculate gradient..."
           let nn = fromIntegral $ step stat
           -- Calculate gradient at current point
           d1s <- lift $ calcGrad (verb params) play nn (dims stat)
           when (verb params) $ info $ "Gra = " ++ show (map gra d1s)
           let -- next point (not adjusted to the limits)
               ds = map (nextPoint nn) d1s
           when (verb params) $ info $ "Nxt = " ++ show (map nxt ds)
           let -- a sequence scaling
               fas osc dim | osc       = (osc, dim)     -- already oscillated
                           | otherwise = ascalStep params nn dim
               (oas, das) = unzip $ zipWith fas (oscs stat) ds
           when (verb params) $ info $ "Alf = " ++ show (map alf das)
           let -- c sequence scaling
               dcs | step stat > mmax params = das
                   | otherwise               = map (cscalStep params nn) das
           when (verb params) $ info $ "Gam = " ++ show (map gam dcs)
           let -- next point
               dns = map (nextStep nn) dcs
           when (verb params) $ info $ "Crt = " ++ show (map crt dns)
           put stat { dims = dns, oscs = oas, step = step stat + 1, grde = grde stat + 1 }
           return False

shiftStep :: Stochastic -> Optim Bool
shiftStep play = do
    stat <- get
    let params = pars stat
    if step stat > nmax params     -- max number of steps
       || (cnxx stat < xstp params && sstp stat >= nstp params)
       then do
           when (verb params) $ do
               if step stat > nmax params
                  then info $ "Termination: nmax steps reached"
                  else info $ "Termination: optimum near enough"
               info $ "Crt = " ++ show (map crt (dims stat))
               info $ "Dims = " ++ show (dims stat)
           return True
       else do
           when (verb params) $ do
               info ""
               info $ "Step " ++ show (step stat) ++ " (shift)"
               info $ "Crt = " ++ show (map crt $ dims stat)
               info $ "calculate gradient..."
           let nn = fromIntegral $ step stat
           -- Calculate gradient at current point
           d1s <- lift $ calcGrad (verb params) play nn (dims stat)
           when (verb params) $ info $ "Gra = " ++ show (map gra d1s)
           let -- next point (not adjusted to the limits)
               ds = map (nextPoint nn) d1s
           when (verb params) $ info $ "Nxt = " ++ show (map nxt ds)
           let -- a sequence shifting
               das | step stat > mmax params = ds
                   | otherwise               = map (bshiftStep params nn) ds
           when (verb params) $ do
               info $ "Bet = " ++ show (map bet das)
               info $ "Vak = " ++ show (map vak das)
           let -- c sequence scaling
               dcs | step stat > mmax params = das
                   | otherwise               = map (cscalStep params nn) das
           when (verb params) $ info $ "Gam = " ++ show (map gam dcs)
           let -- next point, distance in parameter space (for termination)
               xx  = dist1 (map crt dcs) (map nxt dcs)
               dns = map (nextStep nn) dcs
               near | xx < xstp params = sstp stat + 1
                    | otherwise        = 0
           when (verb params) $ do
               info $ "Crt = " ++ show (map crt dns)
               info $ "Dist1 = " ++ show xx
           put stat { dims = dns, step = step stat + 1, grde = grde stat + 1, sstp = near, cnxx = xx }
           return False

info :: String -> Optim ()
info = lift . putStrLn

------------------------
-- Test optimisations --
------------------------

-- 1 dimension, no noise: a simple quadratic funtion with maximum in -10
oneDim :: [Double] -> IO Double
oneDim (x:_) = return $ 1 - 0.1 * (x+10) * (x+10)

maxOneDimExact n = ssSPSA oneDim
                          defSpsaParams { verb = True, nmax = n }
                          [((-50, 100), 50)]

-- Same with noise: noise with a general level
withNoise :: Double -> Stochastic -> [Double] -> IO Double
withNoise noise play xs = do
    y <- play xs
    r <- getStdRandom (randomR (-noise, noise))
    return $ y + r

maxOneDimNoise n = ssSPSA (withNoise 0.5 oneDim)
                          defSpsaParams { verb = True, nmax = n }
                          [((-50, 100), 50)]

-- 2 dimension: negated banana, no noise, (global) maximum at (1, 1)
banana :: [Double] -> IO Double
banana (x:y:_) = return $ negate $ x1 * x1 + 10 * y1 * y1
    where x1 = 1 - x
          x2 = x * x
          y1 = y - x2

maxBanana n = ssSPSA banana
                     defSpsaParams { nmax = n, mmax = n, xstp = 0.001 }
                     [((-10, 10), 5), ((-10, 10), 5)]
