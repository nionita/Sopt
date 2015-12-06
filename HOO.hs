module HOO (
    hooHQ,
    Params(..),
    StartOptions(..),
    hooDefParams,
) where

import Control.Monad

type Stochastic = [Double] -> IO Double
type ORange = (Double, Double)
data StartOptions
    = SOStartWithCheckpoint FilePath [ORange]
    | SOStartNoCheckpoint            [ORange]
    | SORestart             FilePath

-- We represent a region as a list of dimensions with min/max values
-- Here we make integer optimisation, the dimensions have integer type
-- The function which really makes the play is kept with double parameters
-- to remain compatible with the rest of optimisation methods
data HiperQ = HiperQ {
                  dims :: [(Int, Int)],	-- the dimensions of the region
                  func :: Stochastic,	-- function to play once in that region
                  pars :: Params	-- parameters of the algorithm
              }

data Params = Params {
                  verb :: Bool,		-- verbose
                  niu  :: Double	-- smoothness
                  rho  :: Double	-- parameters
              }

hooDefParams :: Params
hooDefParams = Params { verb = False, niu = 10, rho = 0.99 }

hooHQ :: ([Double] -> Double) -> Maybe Params -> StartOptions -> IO [Double]
hooHQ f mp so = do
    case 

divideHyperQ :: HiperQ -> Maybe (HyperQ, HyperQ)
divideHyperQ hq
    | d == 0    = Nothing
    | otherwise = Just (hq1, hq2)
    where lu : lus = dims hq	-- must have at least one dimension
          (i, d) = foldr maxi (0, len lu, 1) lus
          (pre, elm@(l, u) : post) = splitAt i (dims hq)
          med = d `div` 2
          el1 = (l, l + med)
          el2 = (l + med + 1, u)
          hq1 = hq { dims = pre ++ (el1 : post) }
          hq2 = hq { dims = pre ++ (el2 : post) }
          len (l, u) = u - l
          maxi (l, u) (i, ilen, j)
              | jlen > ilen = (j, jlen, j+1)
              | otherwise   = (i, ilen, j+1)
              where jlen = len (l, u)

playHyperQ :: HyperQ -> IO Double
playHyperQ hq = do
    ps <- forM (dims hq) $ \lu -> getStdRandom $ randomR lu
    when (verb $ pars hq) $ putStrLn $ "Crt = " ++ show ps
    func hq $ map fromIntegral ps

instance Region HiperQ where
    divide = divideHyperQ
    play   = playHyperQ

-- A class of regions with 2 methods:
-- one for partitioning a region into 2 disjoint ones
-- one to "play" the bandit given one region, returning the reward
class Region r where
    divide  :: r -> Maybe (r, r)
    rparams :: r -> IO [Double]

-- The HOO tree carries the computation at every step
-- A leaf stores also:
-- - tree level
-- - the current region of the leaf
-- A fork stores also:
-- - tree level
-- - node counter
-- - node estimated mean
-- - node U-value
-- - node B-value
-- - left and right descendents
-- A final node stores:
-- - tree level
-- - node counter
-- - node estimated mean
-- - node U-value
-- - node B-value
-- - the final region (which can't be divided anymore)
data Generic r =>
    HOOTree r = Leaf Int r
              | Fork Int Int Double Double Double (HOOTree r) (HOOTree r)
              | Fin  Int Int Double Double Double r
              deriving Generic

-- Active step: when the tree will play sone arm
astep :: (Region m r, Monad m) => Int -> Double -> Double -> HOOTree r -> m (HOOTree r, Double)
astep n c ro (Leaf h r) = do
    y <- play r
    let k  = 1
        m  = y
        u  = m + sqrt (2 * log n) + c
        b  = u
        h' = h + 1
    case divide r of
        Just (r1, r2) -> return (Fork h k m u b (Leaf h' r1) (Leaf h' r2), y)
        Nothing       -> return (Fin  h k m u b r,                         y)
astep n c ro (Fork h k m u b tl tr) = do
    let (t1, t2) = if tl `better` tr then (tl, tr) else (tr, tl)
        c' = c * ro
    (ta, y) <- astep n c' ro t1
    tb      <- pstep n c' ro t2
    let k' = k + 1		-- new count
        m' = (k * m + y) / k'	-- new mean
        u' = m' + sqrt (2 * log n / k') + c
        b' = minimax u' ta tb	-- min u $ max b1 b2: backward computation
    return (Fork h k' m' u' b' ta tb, y)
astep n c ro (Fin h k m u b r) = do
    y <- play r
    let k' = k + 1		-- new count
        m' = (k * m + y) / k'	-- new mean
        u' = m' + sqrt (2 * log n / k') + c
        b' = minimax u' ta tb	-- min u $ max b1 b2: backward computation
    return (Fin h k' m' u' b' r, y)

-- Passive step: when the tree will not play any arm
pstep :: (Region m r, Monad m) => Int -> Double -> Double -> HOOTree r -> m (HOOTree r)
pstep _ _ _  t = return t
pstep n c ro (Fork h k m u b tl tr) = do
    let c' = c * ro
    ta <- pstep n c' ro tl
    tb <- pstep n c' ro tr
    let u' = m + sqrt (2 * log n / k) + c
        b' = minimax u' ta tb	-- min u $ max b1 b2: backward computation
    return $ Fork h k m u' b' ta tb
pstep n c ro (Fin h k m u b r) = do
    let u' = m + sqrt (2 * log n / k) + c
        b' = u'
    return $ Fin h k m u' b' r

-- Better is the tree with higher b-value
-- Leafs have +infinit there, so they are always better
goodnes :: HOOTree r -> Double
goodnes Leaf {} = maxBound
goodnes (Fork _ _ _ _ b _ _) = b
goodnes (Fin  _ _ _ _ b _)   = b

-- Biased to left tree
better :: HOOTree r -> HOOTree r -> Bool
better l r = goodnes l >= goodnes r

-- Leafs have infinit b-values
minimax :: Double -> HOOTree r -> HOOTree r -> Double
minimax u t1 t2 = min u $ max (goodnes t1) (goodnes t2)

data Generic r =>
    GState = GState {
                 pars :: Params,
                 save :: Maybe FilePath,
                 rfil :: FilePath,
                 step :: !Int,
                 hoot :: HOOTree r
             } deriving Generic

-- Smoothness parameters, a region, a stopping function and an extract function
hoo :: (Region m r, Monad m) => Double -> Double -> r
    -> (Int -> HOOTree r -> m Bool) -> (HOOTree r -> a) -> m a
hoo niu ro a f g = go 1 $ Leaf 0 r
    where go k t = do
              t' <- astep n niu ro t >>= liftM fst
              e <- f k t'
              if e then return (g t') else go (k+1) t'
