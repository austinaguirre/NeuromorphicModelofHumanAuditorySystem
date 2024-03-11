function F = makefilt(M)

F = struct(...
    'b',zeros(M+1,1),...% feedforward filter coefficients
    'a',zeros(M+1,1),...% feedback filter coefficients
    'in',zeros(M,1),...% recent filter inputs 1 (most recent) through M (least recent)
    'out',zeros(M,1),...% recent filter outputs 1 (most recent) through M (least recent)
    'order',M...% filter order
    );

return
%eof