function B0WithoutComma = manipulateB0ValueForSavingPath(B0)

B0String = num2str(B0);
splittedB0String = split(B0String,'.');
B0WithoutComma = [splittedB0String{1} splittedB0String{2}];

end
