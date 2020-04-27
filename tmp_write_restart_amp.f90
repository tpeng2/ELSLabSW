! write restart_amp
program tmp_write_restart_amp
     real amp_save,amp_load
     open(1,file='restart_amp')
     amp_save=0.176410
     amp_load=0.273776
        write(1,*) amp_load,amp_save
     close(1)
end program