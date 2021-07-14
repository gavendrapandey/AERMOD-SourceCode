
f95 -c -O2 -z3 modules.f   
f95 -c -O2 -z3 grsm.f
f95 -c -O2 -z3 aermod.f   
f95 -c -O2 -z3 setup.f    
f95 -c -O2 -z3 coset.f    
f95 -c -O2 -z3 soset.f    
f95 -c -O2 -z3 reset.f    
f95 -c -O2 -z3 meset.f    
f95 -c -O2 -z3 ouset.f    
f95 -c -O2 -z3 inpsum.f   
f95 -c -O2 -z3 metext.f   
f95 -c -O2 -z3 iblval.f   
f95 -c -O2 -z3 siggrid.f  
f95 -c -O2 -z3 tempgrid.f 
f95 -c -O2 -z3 windgrid.f 
f95 -c -O2 -z3 calc1.f    
f95 -c -O2 -z3 calc2.f    
f95 -c -O2 -z3 prise.f    
f95 -c -O2 -z3 prime.f    
f95 -c -O2 -z3 sigmas.f   
f95 -c -O2 -z3 pitarea.f  
f95 -c -O2 -z3 uninam.f 
f95 -c -O2 -z3 output.f   
f95 -c -O2 -z3 evset.f    
f95 -c -O2 -z3 evcalc.f   
f95 -c -O2 -z3 evoutput.f  
f95 -c -O2 -z3 rline.f  
f95 -c -O2 -z3 bline.f 

f95  -o aermod.exe MODULES.obj GRSM.obj AERMOD.obj SETUP.obj COSET.obj SOSET.obj RESET.obj MESET.obj ^
OUSET.obj INPSUM.obj METEXT.obj IBLVAL.obj SIGGRID.obj TEMPGRID.obj WINDGRID.obj CALC1.obj ^
CALC2.obj PRISE.obj PRIME.obj SIGMAS.obj PITAREA.obj UNINAM.obj OUTPUT.obj EVSET.obj EVCALC.obj ^
EVOUTPUT.obj RLINE.obj BLINE.obj

del *.obj
del *.mod