function set_fig(fig, width, height, name)
    % width = 18;
    % height = 7;
    if ~isempty(name)
        set(fig, "Name", name)
    end
    set(fig, "unit", "inches");
    ps = get(fig, "Position");
    set(fig, "Position", [ps(1)/7, ps(2)/7, width, height])
end
