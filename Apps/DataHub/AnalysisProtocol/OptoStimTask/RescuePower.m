function power = RescuePowe(file,p);

    if file == 'Robin061423VIDEOVIEW1.rst'
        power = p;
        power(1:80) = 4.5;
        power(161:160+80)=1.5;
    elseif file == 'Robin062223VIDEOVIEW1.rst'

        power = p;
        power(161:160+80)=6.5;
    elseif file == 'Robin071323VIDEOVIEW1.rst'

        power = p;
        power(481:560)=8.5;

    elseif file == 'Robin072023VIDEOVIEW1.rst'

        power = p;
        power(241:480)=4.5;   

    elseif file == 'Robin072723VIDEOVIEW1.rst'

        power = p;
        power(1601:1760)=1.5;   

 

       

    elseif file =='Adams040324VIDEOVIEW1.rst'
        power = p;
        power(161:240)=1.0;   
    elseif file == 'Adams041624VIDEOVIEW1.rst'
        power = p;
        power(161:240)=0.5;   




    end
end