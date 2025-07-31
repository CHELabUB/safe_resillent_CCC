function handle=figureset(FigNum,varargin)
    % in inches
    defaults={'Height',5.8333,'Width',4.3750};
    options=ScSetOptions(defaults,varargin);
    Height=options.Height;
    Width=options.Width;
    handle=figure(FigNum);
    set(gcf,'units','inches');
    pos_default = get(gcf,'pos');
    pos=pos_default;
    pos(1)=pos(1)-(Width-pos(3))/2;
    pos(2)=pos(2)-(Height-pos(4));
    pos(3)=Width;
    pos(4)=Height;
    % set(gcf,'units','inches','pos',pos,...
    %     'PaperUnits','inches'...
    %     ,'PaperPosition',[0 0 pos(3) pos(4)]);
    set(gcf,'units','inches','pos',pos);
end


function [options,passed_on]=ScSetOptions(defaults,userargs,pass_on)
    %% parses optional arguments of a function and assigns them as fields of
    %% structure options
    % unknown arguments are passed on into cell array passed_on if pass_on is
    % present and non-empty, otherwise an error message is generated
    passed_on={};
    % wrap cell arguments to avoid generating multiple structs
    if isstruct(defaults)
        options=defaults;
    elseif iscell(defaults)
        for i=1:length(defaults)
            if iscell(defaults{i})
                defaults{i}=defaults(i);
            end
        end
        options=struct(defaults{:});
    else
        error('defaults not recognized\n');
    end
    if nargin<3 || isempty(pass_on)
        pass_on=false;
    end
    for i=1:2:length(userargs)
        if isfield(options,userargs{i})
            options.(userargs{i})=userargs{i+1};
        else
            if ~pass_on
                error('option ''%s'' not recognized\n',userargs{i});
            else
                passed_on=[passed_on,{userargs{i},userargs{i+1}}]; %#ok<AGROW>
            end
        end
    end
end
