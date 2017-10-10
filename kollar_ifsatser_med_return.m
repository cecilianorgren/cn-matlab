if 1
    if 1
        a = 0;
        if 1
            reply = input('Do you want more? Y/N [Y]:','s');
            if strcmp(lower(reply),'y')
                disp('Nice!')
            else
                return;
            end
        end
        disp('Return only breaks out of current if-sats.')
        a = 1;
    end
    if a
        disp('Return only breaks out of all if-sats.')
    end
end