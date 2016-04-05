clear
dial = WaitDialog();
try
    i = 1;
    while i < 41
        dial.setMainString(['i: ' mat2str(i)]);
        pause(0.1);
        i = i+1;
        dial.terminate();
    end;
    parfor i=1:40
        k = i;
        dial.terminate();
    end;
    delete(dial);
catch ME
    switch ME.identifier
        case 'MyERR:Terminated'
            warning 'Operation terminated';
        otherwise
            warning 'Rethrowing';
            rethrow(ME);
    end;
end;