f = sprintf('density_save/density%i.gnu',t)
splot f u 1:2:3
t = t+1
pause 0.001
reread 
