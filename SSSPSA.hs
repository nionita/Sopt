module SSSPSA (
    ssSPSA,	-- the optimisation function
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
class Stochastic s where
    play :: s -> [Double] -> IO Double

-- Per dimension needed info
data Dim = Dim {
               low, hig :: !Double,	-- range of dimension
               alf, bet :: !Double,	-- parameters of a(n) sequence
               vak      :: !Double,	-- shift upper bound
               gam      :: !Double,	-- parameter of c(n) sequence
               cmx      :: !Double,	-- maximum for c(n) sequence
               crt      :: !Double,	-- current value
               gra      :: !Double,	-- current gradient value
               nxt      :: !Double,	-- next value
               ash      :: !Int,	-- number of a shifts so far
               csc      :: !Int,	-- number of c sclaes so far
               hil, hih :: Bool	-- hit low/high limit in alpha scaling
           }

data GState = GState {
                 pars :: Params,	-- we include the params in the state record
                 dims :: [Dim],		-- per dimension state
                 oscs :: [Bool],	-- already oscillated dimensions
                 step :: !Int,		-- current step
                 grde :: !Int		-- number of gradient estimations
             }

-- Our monad stack
type Optim a = StateT GState IO a

-- Some parameter of the algorithm (remain constant during one execution):
data Params = Params {
                  phia :: !Double,	-- maximul alpha scaling factor
                  c0   :: !Double,	-- define the maximum gamma: fraction from initial interval
                  gam0 :: !Double,	-- gamma scaling factor
                  gami :: !Double,	-- fraction from initial interval for initial gamma
                  xstp :: !Double,	-- termination when norm1 (xn -xc) < xstp
                  nmax :: !Int,		-- maximum number of steps at all
                  h0   :: !Int,		-- max steps for oscillations (a scaling)
                  ka   :: !Int,		-- maximum number of a(n) shifts per dimension
                  kc   :: !Int,		-- maximum number of c(n) scale ups per dimension
                  gmax :: !Int,		-- max allowed gradients for occilation phase
                  mmax :: !Int		-- max step to which shifts/scalings are allowed
              }

-- Default params
defSpsaParams :: Params
defSpsaParams = Params { phia = 10, c0 = 0.2, gam0 = 2, gami = 20, xstp = 0.5,
                         nmax = 1000, h0 = 4, ka = 50, kc = 50, gmax = 20, mmax = 1000 }

-- Initialisation per dimension
-- Expects: low, high and current (start point)
startDim :: Params -> Double -> Double -> Double -> Dim
startDim par l h c
    = Dim { low = l, hig = h, alf = 1, bet = 0, vak = 10,
            gam = (h - l) / gami par, cmx = c0 par * (h - l), crt = c,
            gra = 0, nxt = 0, ash = 0, csc = 0,
            hil = False, hih = False }

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
                 | otherwise   = (False, alf dim)	-- no scaling when gradient is 0

-- Calculating the alpha' scaling factor for alpha
-- Does not work for gradient 0
{-# INLINE ascal #-}
ascal :: Params -> Double -> Double -> Dim -> (Bool, Double)
ascal par lim n dim
    | f > phia par = (False, alf dim * phia par)	-- limited scaling
    | f < 1        = (False, alf dim)			-- no scaling
    | otherwise    = (True,  alf dim * f)		-- exact scaling
    where f = (lim - crt dim) / ak dim n / gra dim

-- Beta shifting step (per dimension)
-- We expect here that the nxt field is not yet adjusted to the limits
bshiftStep :: Params -> Double -> Dim -> Dim
bshiftStep par n dim
    | ash dim >= ka par                 = dim
    | crt dim <= l && nxt dim > u       = beshi u
    | crt dim >= u && nxt dim < l       = beshi l
    | otherwise                         = dim
    where l = lbk dim n
          u = ubk dim n
          bshift delta = ceiling (alf dim * gra dim / (delta - crt dim) - n - bet dim)
          beshi r = dim { bet = b, vak = v, ash = ash dim + 1 }
              where b' = fromIntegral $ bshift $ r - crt dim
                    (bp, v) | b' > vak dim = (vak dim, vak dim * 2)
                            | otherwise    = (b',      vak dim)
                    b = bet dim + bp

cscalStep :: Params -> Double -> Dim -> Dim
cscalStep par n dim
    | csc dim < kc par &&
      (  crt dim == ubk dim n && nxt dim > crt dim
      || crt dim == lbk dim n && nxt dim < crt dim) = dim { gam = gam dim * gap, csc = csc dim + 1 }
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
calcGrad :: Stochastic s => s -> Double -> [Dim] -> IO [Dim]
calcGrad s n ds = do
    let xs = map crt ds
    dx <- forM ds $ \dim -> do
              d <- randSign
              return $ fromIntegral d * ck dim n
    let xp = zipWith (+) xs dx
        xm = zipWith subtract xs dx
    fp <- play s xp
    fm <- play s xm
    let dgrad dim delta = dim { gra = (fp - fm) / delta }
    return $ zipWith dgrad ds dx

doUntil :: Monad m => m Bool -> m ()
doUntil act = go
    where go = do
              r <- act
              if r then return () else go

ssSPSA :: Stochastic s => s -> Params -> [((Double, Double), Double)] -> IO [Double]
ssSPSA s params dlucs = liftM fst $ runStateT (spsa s) stat
    where stat = GState { pars = params, dims = ds, oscs = os, step = 1, grde = 0 }
          ds = map (\((l, u), c) -> startDim params l u c) dlucs
          os = take (length ds) $ repeat False

spsa :: Stochastic s => s -> Optim [Double]
spsa s = do
    doUntil $ scaleStep s
    doUntil $ shiftStep s
    gets $ map crt . dims

scaleStep :: Stochastic s => s -> Optim Bool
scaleStep s = do
    stat <- get
    let params = pars stat
    if step stat > h0 params || grde stat > gmax params || all id (oscs stat)
       then return True
       else do
           let nn = fromIntegral $ step stat
           -- Calculate gradient at current point
           d1s <- lift $ calcGrad s nn (dims stat)
           let -- next point (not adjusted to the limits)
               ds = map (nextPoint nn) d1s
               -- a sequence scaling
               fas osc dim | osc       = (osc, dim)	-- already oscillated
                           | otherwise = ascalStep params nn dim
               (oas, das) = unzip $ zipWith fas (oscs stat) ds
               -- c sequence scaling
               dcs | step stat > mmax params = das
                   | otherwise               = map (cscalStep params nn) das
               -- next point
               dns = map (nextStep nn) dcs
           put stat { dims = dns, oscs = oas, step = step stat + 1, grde = grde stat + 1 }
           return False

shiftStep :: Stochastic s => s -> Optim Bool
shiftStep s = do
    stat <- get
    let params = pars stat
    if step stat > nmax params	-- max number of steps
       || dist1 (map crt (dims stat)) (map nxt (dims stat)) < xstp params
       then return True
       else do
           let nn = fromIntegral $ step stat
           -- Calculate gradient at current point
           d1s <- lift $ calcGrad s nn (dims stat)
           let -- next point (not adjusted to the limits)
               ds = map (nextPoint nn) d1s
               -- a sequence shifting
               das | step stat > mmax params = ds
                   | otherwise               = map (bshiftStep params nn) ds
               -- next point
               dns = map (nextStep nn) das
           put stat { dims = dns, step = step stat + 1, grde = grde stat + 1 }
           return False
