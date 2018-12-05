% clc, clear, close all;
% mex('-outdir','../toolbox','pq_create.cpp');
% mex('-outdir','../toolbox','pq_delete.cpp');
% mex('-outdir','../toolbox','pq_pop.cpp');
% mex('-outdir','../toolbox','pq_push.cpp');
% mex('-outdir','../toolbox','pq_size.cpp');
% mex('-outdir','../toolbox','pq_top.cpp');
% disp('compile done');


mex pq_create.cpp
mex pq_delete.cpp
mex pq_pop.cpp
mex pq_push.cpp
mex pq_size.cpp
mex pq_top.cpp



