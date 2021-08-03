[![GCC](https://github.com/ricardofrantz/nekStab/actions/workflows/gcc.yml/badge.svg?branch=master)](https://github.com/ricardofrantz/nekStab/actions/workflows/gcc.yml)
[![Intel](https://github.com/ricardofrantz/nekStab/actions/workflows/intel.yml/badge.svg?branch=master)](https://github.com/ricardofrantz/nekStab/actions/workflows/intel.yml)

# Documentation moved to https://ricardofrantz.github.io/nekStabDoc

#Parameter overview:
## DNS
uparam(1) = 0 -> DNS  

## Base-flow

uparam(1) = 1   -> SFD  
uparam(1) = 1.1 -> BoostConv  

## Newton-Krylov

uparam(1) = 2   -> Forced frequency (endTime = period forced)  
uparam(1) = 2.1 -> Natural frequency (endTime = period guess)  

## Eigensolver

uparam(1) = 3.1  -> Direct solver for Steady base flow  
uparam(1) = 3.11 -> Direct solver for Periodic base flow  
uparam(1) = 3.2  -> Adjoint solver for Steady base flow  
uparam(1) = 3.21 -> Adjoint solver for Periodic base flow  
uparam(1) = 3.3  -> Transient Growth  

## Postprocessing
uparam(1) = 4          -> Perturbation kinetic energy computation  
uparam(1) = 4.1        -> Wavemaker computation  
uparam(1) = 4.2        -> Base flow sensitiviy computation  
uparam(1) = 4.31/4.32  -> Forcing sensitiviy computation  
