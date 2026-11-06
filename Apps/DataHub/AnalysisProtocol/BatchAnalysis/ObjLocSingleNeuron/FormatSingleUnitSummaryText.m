function L = FormatSingleUnitSummaryText(PD, SA, SB, Cmp, names)
% FormatSingleUnitSummaryText  Plain-text observations + area comparison,
% built from the computed results (for the console / a .txt for slides).
%
% Inputs : PD (ComputePresentationData), SA, SB (ComputeSingleUnitSummary),
%          Cmp (ComputeAreaComparison), names {'FEF','SC'}
% Output : L — cellstr, one line per entry
%
% Xuefei Yu, 2026-11-05

S = {SA, SB}; L = {};
L{end+1} = '=== Single-unit summary: object (stable) vs location (oneDR) value ===';
L{end+1} = sprintf('Neurons: %s %d, %s %d (same selection as Organize: FR>1, visual, >=10 trials/condition, opto)', ...
    names{1}, SA.n, names{2}, SB.n);
L{end+1} = '';
L{end+1} = sprintf('1. %s', PD.Fig1.Headline);
for j = 1:2
    L{end+1} = sprintf('   %s: object %.0f%%, location %.0f%% of neurons (chance 5%%)', names{j}, 100*PD.Fig1.Frac(1,j), 100*PD.Fig1.Frac(2,j)); %#ok<AGROW>
end
L{end+1} = sprintf('   Area difference: object %s; location %s', SlidePString(PD.Fig1.pArea(1), names), SlidePString(PD.Fig1.pArea(2), names));
L{end+1} = sprintf('2. %s', PD.Fig2.Headline);
for j = 1:2
    L{end+1} = sprintf('   %s: pre-target difference in object task %.0f%%, location task %.0f%%; same sign as response %.0f%%', ...
        names{j}, 100*PD.Fig2.Frac(1,j), 100*PD.Fig2.Frac(2,j), 100*S{j}.BaseRespAgree.Fraction(1)); %#ok<AGROW>
end
L{end+1} = sprintf('   Area difference (location task): %s', SlidePString(PD.Fig2.pArea(2), names));
L{end+1} = sprintf('3. %s', PD.Fig3.Headline);
for j = 1:2
    L{end+1} = sprintf('   %s: both-coding %d vs %.1f expected, OR %.2f [%.2f %.2f]', names{j}, PD.Fig3.Obs(j), PD.Fig3.Exp(j), ...
        PD.Fig3.Est(j), PD.Fig3.Lo(j), PD.Fig3.Hi(j)); %#ok<AGROW>
end
L{end+1} = sprintf('4. %s', PD.Fig5.Headline);
for j = 1:2
    L{end+1} = sprintf('   %s: observed r %.2f [%.2f %.2f], reliability obj %.2f / loc %.2f, noise-corrected r %.2f [%.2f %.2f]', ...
        names{j}, PD.Fig5.rObs(j), PD.Fig5.rObsCI(j,:), PD.Fig5.relObj(j), PD.Fig5.relLoc(j), PD.Fig5.rDis(j), PD.Fig5.rDisCI(j,:)); %#ok<AGROW>
end
L{end+1} = sprintf('   Area difference (noise-corrected r): %s', SlidePString(PD.Fig5.pArea, names));
for j = 1:2
    c = S{j}.Concord;
    L{end+1} = sprintf('   %s same preferred value: raw %.0f%%, corrected %.0f%%', names{j}, 100*c.Raw.Fraction, 100*c.Corrected.Fraction); %#ok<AGROW>
end
L{end+1} = '';
sig = Cmp(Cmp.pHolm < 0.05, :);
if isempty(sig)
    L{end+1} = sprintf('FEF vs SC: no measure differs after Holm correction (%d tests).', height(Cmp));
else
    L{end+1} = sprintf('FEF vs SC: %d of %d measures differ after Holm correction:', height(sig), height(Cmp));
    for i = 1:height(sig)
        L{end+1} = sprintf('   %s: %s %.3f vs %s %.3f (%s, p_Holm = %.3g)', sig.Measure{i}, names{1}, sig.(names{1})(i), ...
            names{2}, sig.(names{2})(i), sig.Test{i}, sig.pHolm(i)); %#ok<AGROW>
    end
end
end
