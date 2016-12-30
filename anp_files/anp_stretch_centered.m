% stretches a given interval with a given factor about its center
function result = anp_stretch_centered(in, factor)
    size_in = in(2) - in(1);
    result = [(in(1)-(factor-1)*size_in/2),(in(2)+(factor-1)*size_in/2)];
end
