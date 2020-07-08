function CW=Increase(CW,CWmax)
%增大退避窗口大小
CW=2*CW;
if CW>=CWmax
    CW=CWmax;
end
end

