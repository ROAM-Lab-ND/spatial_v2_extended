function out = outputSelect(varargin)
    num = varargin{1};
    f = varargin{2};
    eval_str = 'f(';
    for i = 3:(nargin)
        eval_str = [eval_str 'varargin{' num2str(i) '},'];
    end
    eval_str = [eval_str(1:end-1) ')'];

    assert(num <= 8,'requires num <= 8');
    
    if num ==1
        out = eval(eval_str);
    elseif num==2
        [~, out] = eval(eval_str);
    elseif num==3
        [~,~,out] = eval(eval_str);
    elseif num==4
        [~,~,~,out] = eval(eval_str);
    elseif num==5
        [~,~,~,~,out] = eval(eval_str);
    elseif num==6
        [~,~,~,~,~,out] = eval(eval_str);
    elseif num==7
        [~,~,~,~,~,~,out] = eval(eval_str);
    elseif num==8
        [~,~,~,~,~,~,~,out] = eval(eval_str);
    end
end