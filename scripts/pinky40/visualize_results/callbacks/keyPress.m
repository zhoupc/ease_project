
function keyPress(~, e)
switch lower(e.Key)
    case 'n' % check enxt cell
        evalin('base', 'cell_id = cell_id+1; show_neuron;');
    case 'b' % check enxt cell
        evalin('base', 'cell_id = cell_id-1; show_neuron;');
end