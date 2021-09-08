function B0WithoutComma = manipulateB0ValueForSavingPath(B0)

B0String = num2str(B0);
try
    splittedB0String = split(B0String,'.');
    B0WithoutComma = [splittedB0String{1} splittedB0String{2}];
catch
    B0WithoutComma = B0String;
    disp('B0 has no comma.')
end


end
