%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script to create compiled code for stand-alone matlab run
% Still in progress, but the basic idea is to create a working 
% directory on scratch where the code required is moved to, and
% and use it for compilation and running the new code.
% I clumsily address the problem of specifying function dependencies by
% making sure the main function paths are removed - the only path is local
% this is done by a separate function that is substituted before compilation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 do_compile = 0;	% 0 compilation done outside matlab, require call
			% 1 compilation done by this function
 curDir = pwd;
 scratchDir = '/u/scratch2/scratch1/d/danieleb/';
 baseDir = 'test1/';
 workDir = [scratchDir baseDir]; 

 % Main script to compile
 mainFun = 'boats2d_suite_biomass3_par.m';

 % Compilation options
 % -m : specifies C language translation during compilation
 % -R : specifies runtime options
%compOpt = ['-m -R -singleCompThread -R -nosplash -R -nojvm']; 
 compOpt = ['-m -R -nosplash']; 

 % Directory for compilation
 funName = mainFun(1:end-2);
 compDir = [workDir funName];

 % Create directories for compilation
 if ~(exist(workDir)==7)
    mkdir(workDir)
 end
 if ~(exist(compDir)==7)
    mkdir(compDir)
 end

 % Make sure local paths are present
 add_local_paths

 % Find file dependencies
 fList = matlab.codetools.requiredFilesAndProducts(mainFun);

 % Copies all file dependencies to compilation directory
 nList = length(fList);
 for indf=1:nList
   copyfile(fList{indf},compDir);
 end

 cd(compDir)

 % Removes (or substitute) paths definitions  
 cstring = ['!echo "path = ' '''' './' '''' ';" > add_local_paths.m'];
 eval(cstring);
 
 % Updates function list
 tDir = dir;
 tList = {};
 for indf=1:length(tDir)
    if strfind(fliplr(tDir(indf).name),'m.')==1
       tList = [tList tDir(indf).name];
    end
 end
 tList = setdiff(tList,mainFun);
 
 % Prepare function call for compiler - can add options here
 sComp = ['mcc ' compOpt ' '  mainFun];
 for indf=1:length(tList)
    sComp = [sComp ' ' tList{indf}];
 end

 % Evaluate/display compiling command
 switch do_compile
 case 0
    disp(['To create and run compiled ' mainFun ' :']);
    disp(['cd ' compDir]);
    disp(['module load matlab']);
    disp(sComp)
    disp(['./' funName]);
 case 1
    eval(sComp)
    disp(['To run compiled ' mainFun ' :']);
    disp(['cd ' compDir]);
    disp(['module load matlab']);
    disp(['./' funName]);
 otherwise
    error(['do_compile option not valid'])
 end

 % Back to current directory
 cd(curDir); 

