{-# LANGUAGE DeriveGeneric #-}
{-# LANGUAGE MultiParamTypeClasses #-}
{-# LANGUAGE FunctionalDependencies #-}

module Optimisation.Stochastic.Common (
    Stochastic,
    ORange,
    StartOptions(..),
    OptimState(..),
    Dimension(..),
    Optim,
    optimise,
    info,
    checkEnd,
    checkRunning,
    calcNGrad
  ) where

import Control.Monad
import Control.Monad.IO.Class
import Control.Monad.Trans.Class
import Control.Monad.Trans.State.Strict
import Control.Monad.Trans.Reader
import qualified Data.ByteString.Lazy as BL
import Data.Default
import Data.Serialize (Serialize, encodeLazy, decodeLazy)
import GHC.Generics
import System.Directory
import System.FilePath
import System.Random

-- We need to let some external process to get a sample of the function
-- at some position, and this will happen in IO
type Stochastic = [Double] -> IO Double
type ORange = (Double, Double)
data StartOptions
    = SOStartWithCheckpoint FilePath [ORange]
    | SOStartNoCheckpoint            [ORange]
    | SORestart             FilePath

-- A class for a state of an optimisation
-- It has an associated parametr type and some methods to set/query
-- some important information, like current point, save file
class (Default s, Default p, Serialize s, Serialize p)
    => OptimState s p | p -> s, s -> p where
    iniCurr :: [ORange] -> s -> s
    getCurr :: s -> [Double]
    setSave :: FilePath -> s -> s
    getSave :: s -> Maybe FilePath
    setRun  :: FilePath -> s -> s
    getRun  :: s -> Maybe FilePath

-- Class to define some operations per dimension, used
-- to calcuate the gradients
class Dimension d where
    dimCrt :: d -> Double       -- current value of dimension
    dimDlt :: d -> Double       -- delta for gradient
    setGrd :: Double -> d -> d  -- set the current grandient value

-- Our monad stack with parameters and state
type Optim s p a = ReaderT p (StateT s IO) a

-- Some utilities:
doUntil :: OptimState s p => Optim s p Bool -> Optim s p ()
doUntil act = go
    where go = do
              r <- act
              if r then
                 return ()
                 else do
                     checkPoint
                     go

runFileNum :: IO Int
runFileNum = getStdRandom (randomR (1000, 9999))

optimise :: OptimState s p
         => Stochastic -> (Stochastic -> Optim s p Bool) -> Bool -> Maybe p -> StartOptions
         -> IO [Double]
optimise play ostep rf mparams staopts = do
    (stat, params)
        <- case staopts of
                SOStartWithCheckpoint file dlucs -> do
                    let params = maybe def id mparams
                        stat   = setSave file $ iniCurr dlucs def
                    return (stat, params)
                SOStartNoCheckpoint        dlucs -> do
                    let params = maybe def id mparams
                        stat   = iniCurr dlucs def
                    return (stat, params)
                SORestart             file       -> do
                    (stat', pars') <- restorePoint file
                    let params = maybe pars' id mparams  -- we can override params
                        stat   = setSave file stat'
                    return (stat, params)
    -- If requested, write an empty file with a random number in name to be used for stopping
    -- an optimisation run as soon as possible, but in a smooth way
    -- If the file will be deleted, the optimisation will stop after the current step ends
    -- If you have no checkpoint file, you will not be able to continue from there!
    stat' <- if rf
                then do
                    r <- runFileNum
                    let rfile = "running-" ++ show r
                    BL.writeFile rfile BL.empty
                    return $ setRun rfile stat
                else return stat
    (ds, _) <- runStateT (runReaderT (optim play ostep) params) stat
    return ds

optim :: OptimState s p => Stochastic -> (Stochastic -> Optim s p Bool) -> Optim s p [Double]
optim play ostep = do
    doUntil $ ostep play
    lift $ gets getCurr

info :: OptimState s p => String -> Optim s p ()
info = liftIO . putStrLn

checkEnd :: OptimState s p => Bool -> Bool -> String -> Optim s p Bool -> Optim s p Bool
checkEnd verb cond mes act
    | cond      = do
        when verb $ info mes
        return True
    | otherwise = act

checkPoint :: OptimState s p => Optim s p ()
checkPoint = do
    s <- lift get
    case getSave s of
        Nothing -> return ()
        Just cp -> liftIO $ do
            let cpn = addExtension cp "new"
            BL.writeFile cpn $ encodeLazy s
            renameFile cpn cp

restorePoint :: OptimState s p => FilePath -> IO (s, p)
restorePoint file = do
    bs <- BL.readFile file
    case decodeLazy bs of
        Left estr -> error $ "Decode: " ++ estr
        Right sta -> return sta

checkRunning :: OptimState s p => Optim s p Bool -> Optim s p Bool
checkRunning act = do
    mrf <- lift $ gets getRun
    case mrf of
        Nothing -> act
        Just rf -> do
            rfe <- liftIO $ doesFileExist rf
            if rfe then act else do
               liftIO $ putStrLn $ "Running file " ++ rf ++ " not found, stopping"
               -- put s { rfdl = False }
               return True

-- Uniformly random +/-1
randSign :: IO Int
randSign = do
    r <- getStdRandom (randomR (1,2))
    case r of
        1 -> return r
        2 -> return (-1)

-- Calculate k gradients
-- We are stochastic, so taking a mean of a few "gradients"
-- could be beneficial, that's why k is
calcNGrad :: Dimension d => Bool -> Stochastic -> Int -> [d] -> IO [d]
calcNGrad verb play k ds = do
    let xs = map dimCrt ds
    dxs <- sequence $ replicate k $ forM ds $ \dim -> do
               d <- randSign
               return $ fromIntegral d * dimDlt dim
    let xps = map (zipWith (+)             xs) dxs
        xms = map (zipWith (flip subtract) xs) dxs
    grs <- forM (zip xps xms) $ \(xp, xm) -> do
        fp <- play xp
        when verb $ putStrLn $ "Func + at " ++ show xp ++ ": " ++ show fp
        fm <- play xm
        when verb $ putStrLn $ "Func - at " ++ show xm ++ ": " ++ show fm
        let df = fp - fm
            dx = zipWith subtract xm xp
        return $ map (\dx -> df / dx) dx
    let invk = 1 / (fromIntegral k)
        gs   = map (*invk) $ foldr1 (zipWith (+)) grs
        dgs  = zipWith (\dim g -> setGrd g dim) ds gs
    return dgs
