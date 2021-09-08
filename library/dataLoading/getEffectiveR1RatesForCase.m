function [effectiveLipidR1Rates,effectiveWaterR1Rates] = ...
    getEffectiveR1RatesForCase(caseForEffectiveR1Rates,lipidData,waterData)

caseForEffectiveR1Rates = lower(caseForEffectiveR1Rates);
switch caseForEffectiveR1Rates
    case 'mean'
        effectiveLipidR1Rates = lipidData.effectiveRelaxationRatesMean;
        effectiveWaterR1Rates = waterData.effectiveRelaxationRatesMean;
    case 'median'
        effectiveLipidR1Rates = lipidData.effectiveRelaxationRatesMedian;
        effectiveWaterR1Rates = waterData.effectiveRelaxationRatesMedian;
    otherwise
        error(['Unknown case for effective R1 rates.' ...
            'Take a look at the configuration file.'])
end
end