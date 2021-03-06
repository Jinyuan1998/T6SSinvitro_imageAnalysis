function T = TimeLapseCellKilling(folderName, useI, printH)
%%%%% **********************************************************************************%%%%% 
%%%%%  image analysis of time lapse using 63x oil, fluorescent cells.                   %%%%% 
%%%%%  bpass function is from the track package by David Grier, John                    %%%%%
%%%%%  Crocker, and Eric Weeks.                                                         %%%%%
%%%%%  (https://site.physics.georgetown.edu/matlab/code.html)                           %%%%%
%%%%% **********************************************************************************%%%%%

%%%%% **********************************************************************************%%%%% 
%%%%%  This matlab script is used to track Aeromonas veronii cell killed                %%%%% 
%%%%%  by the type VI secretion system of Vibrio cholerae in vitro.                     %%%%%
%%%%%  https://www.biorxiv.org/content/early/2017/11/29/226472.full.pdf+html            %%%%%
%%%%% **********************************************************************************%%%%%

%%% to do: the predator track still needs work due to weak fluoresent

% folderName: the directory containing time lapse image
% useI: the index of images used for analysis
% printH: handle controls whether to print image or not; default is 0=not print; 1=print


cd(folderName)
f = dir('Experiment*');
F = f(find(vertcat(f.isdir)));
folder = F.name;
% image analysis of the time lapse
% output the number of prey/predator cells in each frame of the time lapse
% predator: constitutive vibrio killer (green) + prey: Aeromonas veronii (orange)
useimage = useI;
%% 1. check to see if the subfolders are made
if ~ exist('outputPD', 'dir')
    mkdir outputPD
end
if ~ exist('OutputPrey', 'dir')
    mkdir OutputPrey
end
if ~ exist('cropPd', 'dir')
    mkdir cropPd
end
if ~ exist('cropPrey', 'dir')
    mkdir cropPrey
end
if ~ exist('prcellplot', 'dir')
    mkdir prcellplot
end

preyimagefiles = dir([folder '/*Cy3*.tif']);
pdtimagefiles = dir([folder '/*EGFP*.tif']);
preyimagefiles = sort({preyimagefiles.name});
pdtimagefiles = sort({pdtimagefiles.name});

printFigure = printH;

%% 2. use tracking package by Daniel Blair and Eric Dufresne
% (http://site.physics.georgetown.edu/matlab/code.html)
% only analyze predator for now
for i=useimage
    pr = double(imread([folder '/' preyimagefiles{i}]));    
    prb = bpass(pr, 2, 5);
    iIdx = i;
    if iIdx > 1  % the view of microscope drifted during time lapse; now align those figures using the first as reference
        prfix = double(imread([folder '/' preyimagefiles{i-1}]));
        prfix = bpass(prfix, 1, 5);
        
        [optimizer, metric] = imregconfig('multimodal');        
        tform = imregtform(prb,prfix,'translation',optimizer,metric);
        movingRegistered = imwarp(prb,tform,'OutputView',imref2d(size(prfix)));
        if printFigure == 1
            figure; imshowpair(imadjust(prfix),imadjust(movingRegistered),'Scaling','joint')
        end
        % get the coordinates that transforms every time
        X(i) = tform.T(3,1);
        Y(i) = tform.T(3,2);        
    end
end
X(X~=0) 
Y(Y ~= 0)

%% 3. crop drifted images
X1 = sum(X(X>=0));
X2 = sum(X(X<0));
Y1 = sum(Y(Y>=0));
Y2 = sum(Y(Y < 0));
[H, W] = size(imread([folder '/' preyimagefiles{1}]));
Wc = round(W - abs(X1) - abs(X2))-2;
Hc = round(H - abs(Y1) - abs(Y2))-2;

for i=useimage
    %read raw images
    pr = imread([folder '/' preyimagefiles{i}]);
    pd = imread([folder '/' pdtimagefiles{i}]);
    
    xs = max([sum(X(i+1:end)), 0, -sum(X(1:i))]);
    ys = max([sum(Y(i+1:end)), 0, -sum(Y(1:i))]);
    pr_c = imcrop(pr,[round(xs)+1 round(ys)+1 Wc Hc]);
    pd_c = imcrop(pd,[round(xs)+1 round(ys)+1 Wc Hc]);

    p1 = preyimagefiles{i};
    imwrite(pr_c, ['cropPrey/'  p1(1:end-4) '_crop.tif']);
    p2 = pdtimagefiles{i};
    imwrite(pd_c, ['cropPd/' p2(1:end-4) '_crop.tif']);
       
    imwrite(imadjust(pr_c), ['outputPrey/Crop'  p1(1:end-4) '.tif']);
    imwrite(imadjust(pd_c), ['outputPd/Crop'  p2(1:end-4) '.tif']);

end

%% 4.find particles and count
close all
ctpreyimagefiles = dir('cropPrey/*_crop.tif');
ctpdtimagefiles = dir('cropPd/*_crop.tif');
ctpreyimagefiles = sort({ctpreyimagefiles.name});
ctpdtimagefiles = sort({ctpdtimagefiles.name});
for i=1:length(useimage)
    pr = imread(['cropPrey/' ctpreyimagefiles{i}]);    
    prb = bpass(pr, 2, 5, median(median(pr)) / ((200-useimage(i))/3));
            
    % use watershed algorith %%%%%%%%%%  
    prbw = im2bw(imadjust(prb));
    prbw = bwareaopen(prbw, round(60/(1+i/10)));
    L = watershed(prbw);
    L(~prbw) = 0;
    B = uint16(L) .* pr;
    bb = bwconncomp(B);
    Lbb = labelmatrix(bb);
    PRrgb = label2rgb(Lbb, 'spring', 'c', 'shuffle'); 
    if printFigure == 1
        figure; imagesc(imadjust(pr))
        figure; imshow(PRrgb)
    end
    prkw = pkfnd(B, 0.01, 11);
    prcnt = cntrd(double(Lbb), prkw, 4);
    warning off
    prpos_list = prcnt(:,1:2);
    prcountws(i) = size(prpos_list, 1);
    clear p
    
    %%%%%%%%%%%%%%%%%%%%%%%%% analyze predator (not working well for now due to the weak fluorescence) %%%%%%%%%%%%%%%%
    pd = imread(['cropPd/' ctpdtimagefiles{i}]);
    pdb = bpass(pd, 1, 4, median(median(pd))/(200-useimage(i)))-prb;
    pdb(pdb<0) = 0;
    pdbw = im2bw(imadjust(pdb));
    pdbw = bwareaopen(pdbw, 40);
    
    pdL = watershed(pdbw, 4);
    pdL(~pdbw) = 0;
    C = uint16(pdL) .* pd;
    cc = bwconncomp(C);
    Lcc = labelmatrix(cc);
    PDrgb = label2rgb(Lcc,'spring','k','shuffle');

    pdkw = pkfnd(pdbw, 0.01, 10);
    pdcnt = cntrd(pdbw,pdkw,11);
    pdpos_list = pdcnt(:,1:2);
    pdcountws(i) = size(pdpos_list, 1);
    if printFigure == 1
        figure; imagesc(imadjust(pd))
        figure; imshow(PDrgb)
        figure, imshowpair(prbw, pdbw, 'Scaling','joint'); 
%         hold on; plot(pdpos_list(:,1), pdpos_list(:,2), 'yx');
        hold on; plot(prkw(:,1), prkw(:,2), 'r*')
        hold off
    end
%     figure
%     imshow(imadjust(pr));
%     hold on; plot(prpos_list(:,1), prpos_list(:,2), 'rx', 'MarkerSize', 8);
%     hold off
%     print(['prcellplot/' ctpreyimagefiles{i}], '-dtiff')
%     close
        
end

A = [useimage' prcountws' pdcountws'];
T = array2table(A, 'variableNames', {'imageN' 'prey_N' 'predator_N'});
