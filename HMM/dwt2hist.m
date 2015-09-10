function dwt2hist(BIN,LEN)
% USAGE: dwt2hist(BIN,LEN)
% BIN is the time resolution in msec
% LEN is the length of the traces in frames

%   Copyright 2007-2015 Cornell University All Rights Reserved.


if nargin == 0
    BIN=10; %ms
    LEN=1500; %frames
end

LEN=LEN*BIN;

zero2one=[];
zero2two=[];
zero2three=[];
one2zero=[];
one2two=[];
one2three=[];
two2zero=[];
two2one=[];
two2three=[];
three2zero=[];
three2one=[];
three2two=[];

trans_no=[];

time=0:BIN:LEN;

[filename, pathname]=uigetfile('*.dwt','Select QuB dwell times file.');
if filename==0
    return
end
file=strcat(pathname,filename);
fid=fopen(file,'r');

while 1
    dwells=[];

    line=fgetl(fid);
    if ~ischar(line)
        break;
    end

    if line(1)=='S'
        line=fgetl(fid);
    end

    %---Make a list of dwells from one trace at a time
    while line(1)~='S'
        dwells=[dwells; strread(line)];
        line=fgetl(fid);

        if ~ischar(line)
            break;
        end
    end

    seq=dwells(:,1);
    ndwells=numel(seq);

    trans_no=[trans_no; ndwells];

    for m=2:ndwells
        switch seq(m-1)
            case 0
                switch seq(m)
                    case 1
                        zero2one=[zero2one; dwells(m-1,2)];
                    case 2
                        zero2two=[zero2two; dwells(m-1,2)];
                    case 3
                        zero2three=[zero2three; dwells(m-1,2)];
                end
            case 1
                switch seq(m)
                    case 0
                        one2zero=[one2zero; dwells(m-1,2)];
                    case 2
                        one2two=[one2two; dwells(m-1,2)];
                    case 3
                        one2three=[one2three; dwells(m-1,2)];
                end
            case 2
                switch seq(m)
                    case 0
                        two2zero=[two2zero; dwells(m-1,2)];
                    case 1
                        two2one=[two2one; dwells(m-1,2)];
                    case 3
                        two2three=[two2three; dwells(m-1,2)];
                end
            case 3
                switch seq(m)
                    case 0
                        three2zero=[three2zero; dwells(m-1,2)];
                    case 1
                        three2one=[three2one; dwells(m-1,2)];
                    case 2
                        three2two=[three2two; dwells(m-1,2)];
                end
        end

    end

end

fclose(fid);

h_zero2one=hist(zero2one,time);
h_zero2two=hist(zero2two,time);
h_zero2three=hist(zero2three,time);
h_one2zero=hist(one2zero,time);
h_one2two=hist(one2two,time);
h_one2three=hist(one2three,time);
h_two2zero=hist(two2zero,time);
h_two2one=hist(two2one,time);
h_two2three=hist(two2three,time);
h_three2zero=hist(three2zero,time);
h_three2one=hist(three2one,time);
h_three2two=hist(three2two,time);

c_zero2one=cumsum(h_zero2one)';
c_zero2two=cumsum(h_zero2two)';
c_zero2three=cumsum(h_zero2three)';
c_one2zero=cumsum(h_one2zero)';
c_one2two=cumsum(h_one2two)';
c_one2three=cumsum(h_one2three)';
c_two2zero=cumsum(h_two2zero)';
c_two2one=cumsum(h_two2one)';
c_two2three=cumsum(h_two2three)';
c_three2zero=cumsum(h_three2zero)';
c_three2one=cumsum(h_three2one)';
c_three2two=cumsum(h_three2two)';

s_zero2one=(max(c_zero2one)-c_zero2one)/max(c_zero2one);
s_zero2two=(max(c_zero2two)-c_zero2two)/max(c_zero2two);
s_zero2three=(max(c_zero2three)-c_zero2three)/max(c_zero2three);
s_one2zero=(max(c_one2zero)-c_one2zero)/max(c_one2zero);
s_one2two=(max(c_one2two)-c_one2two)/max(c_one2two);
s_one2three=(max(c_one2three)-c_one2three)/max(c_one2three);
s_two2zero=(max(c_two2zero)-c_two2zero)/max(c_two2zero);
s_two2one=(max(c_two2one)-c_two2one)/max(c_two2one);
s_two2three=(max(c_two2three)-c_two2three)/max(c_two2three);
s_three2zero=(max(c_three2zero)-c_three2zero)/max(c_three2zero);
s_three2one=(max(c_three2one)-c_three2one)/max(c_three2one);
s_three2two=(max(c_three2two)-c_three2two)/max(c_three2two);

outfile=strrep(file,'.dwt','_surv.txt');

fid2=fopen(outfile,'w');
fprintf(fid2,'%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n',...
    'time','zero2one','zero2two','zero2three',...
    'one2zero','one2two','one2three',...
    'two2zero','two2one','two2three',...
    'three2zero','three2one','three2two');

n_zero2one=numel(s_zero2one);
n_zero2two=numel(s_zero2two);
n_zero2three=numel(s_zero2three);
n_one2zero=numel(s_one2zero);
n_one2two=numel(s_one2two);
n_one2three=numel(s_one2three);
n_two2zero=numel(s_two2zero);
n_two2one=numel(s_two2one);
n_two2three=numel(s_two2three);
n_three2zero=numel(s_three2zero);
n_three2one=numel(s_three2one);
n_three2two=numel(s_three2two);

rows=max([n_zero2one n_zero2two n_zero2three...
    n_one2zero n_one2two n_one2three...
    n_two2zero n_two2one n_two2three...
    n_three2zero n_three2one n_three2two]);

%Write each lifetime array to a different column, padding with '-'
for k=1:rows

    fprintf(fid2,'%g',time(k)/1000);
    fprintf(fid2,'\t');
    
    if k<=n_zero2one, fprintf(fid2,'%g\t',s_zero2one(k,1));
    else fprintf(fid2,'%s\t','--'); end

    if k<=n_zero2two, fprintf(fid2,'%g\t',s_zero2two(k,1));
    else fprintf(fid2,'%s\t','--'); end

    if k<=n_zero2three, fprintf(fid2,'%g\t',s_zero2three(k,1));
    else fprintf(fid2,'%s\t','--'); end

    if k<=n_one2zero, fprintf(fid2,'%g\t',s_one2zero(k,1));
    else fprintf(fid2,'%s\t','--'); end

    if k<=n_one2two, fprintf(fid2,'%g\t',s_one2two(k,1));
    else fprintf(fid2,'%s\t','--'); end

    if k<=n_one2three, fprintf(fid2,'%g\t',s_one2three(k,1));
    else fprintf(fid2,'%s\t','--'); end

    if k<=n_two2zero, fprintf(fid2,'%g\t',s_two2zero(k,1));
    else fprintf(fid2,'%s\t','--'); end

    if k<=n_two2one, fprintf(fid2,'%g\t',s_two2one(k,1));
    else fprintf(fid2,'%s\t','--'); end

    if k<=n_two2three, fprintf(fid2,'%g\t',s_two2three(k,1));
    else fprintf(fid2,'%s\t','--'); end

    if k<=n_three2zero, fprintf(fid2,'%g\t',s_three2zero(k,1));
    else fprintf(fid2,'%s\t','--'); end

    if k<=n_three2one, fprintf(fid2,'%g\t',s_three2one(k,1));
    else fprintf(fid2,'%s\t','--'); end

    if k<=n_three2two, fprintf(fid2,'%g\n',s_three2two(k,1));
    else fprintf(fid2,'%s\n','--'); end

end
fclose(fid2);
