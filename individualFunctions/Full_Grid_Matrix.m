function out=Full_Grid_Matrix(in)

%create a full grid matrix of input variables
%used to create the full grid of all txffe cases
%example:
%Full_Grid_Matrix({ [1 2] [100 200] })
%out =
%     1   100
%     1   200
%     2   100
%     2   200
%
%input can also be mixed between numeric and cell of char
%example:
%Full_Grid_Matrix({ [1 2] {'A' 'B'} })
%out =
%    {[1]}    {'A'}
%    {[1]}    {'B'}
%    {[2]}    {'A'}
%    {[2]}    {'B'}

if ~iscell(in)
    error('input must be cell array of individual sweep variables');
end

num_columns=length(in);
num_cases=prod(cellfun('length',in));

cell_output=0;
cell_indices=cellfun(@(x) iscell(x),in);
if any(cell_indices)
    cell_output=1;
end
if cell_output
    for k=find(~cell_indices)
        in{k}=num2cell(in{k});
    end
end

if cell_output
    out=cell(num_cases,num_columns);
else
    out=zeros(num_cases,num_columns);
end

%num_repetitions controls how many times each element of the column
%repeats.  The first column is always just a copy of itself since every
%case will vary.
num_repetitions=1;
for k=num_columns:-1:1
    this_column=in{k}(:);
    %copy the column into a matrix to create the repetitions needed
    B=repmat(this_column,[1 num_repetitions]);
    %reshape into single column (actual repetitions)
    C=reshape(B',[numel(B) 1]);
    %repeat the single column to build the entire length required
    num_repeats=num_cases/length(C);
    D=repmat(C,[num_repeats 1]);
    out(:,k)=D;
    %determine how many repetitions the next column needs
    num_repetitions=num_repetitions*length(this_column);
end