clear all
% 对于40和50 都挺不错
%% 基本参数设定
KyN=180;
KzN=40;

NE=20;
%centerRatio=linear_map(NE);
centerRatio=0.3;
%% 得到小椭圆
a=double(KyN/2.0);
b=double(KzN/2.0);
ystart=int32(-a);
yend=int32(a);
ycenter=int32(a);

zstart=int32(-b);
zend=int32(b);
zcenter=int32(b);

a2_2=a*a;
b2_2=b*b;

lambda1 = 4;
lambda2 =12;
a1_2 = a2_2 / lambda1;
b1_2 = b2_2 / lambda2;
% 确保小椭圈的面积为大椭圆的centerRatio左右
while (abs(sqrt(a1_2) * sqrt(b1_2)  / ( a  *  b ) -centerRatio) > 0.05)
    if (sqrt(a1_2) * sqrt(b1_2) / (a * b) > centerRatio)

        a1_2 = a1_2 * 0.95;
        b1_2 = b1_2 * 0.95;

    else

        a1_2 = a1_2 * 1.05;
        b1_2 = b1_2 * 1.05;
    end
end
TrajIn=[];
TrajOut=[];
TrajPeri=[];
for ky = double(ystart: yend-1)
    for kz = double(zstart:zend-1)
        %idx = (ky + ycenter) * KzN + (kz + zcenter);
        % 椭圆方程
        tmp1 = ky * ky / a1_2 + kz * kz / b1_2;
        tmp2 = ky * ky / a2_2 + kz * kz / b2_2;
        % 椭圆面积S=π(圆周率)×a×b
        if (tmp1 < 1)% 得到椭圆内部从中心开始的25%的点
            TrajIn=cat(1,TrajIn,[ky,kz]); % The central K-space sector(inner circle)
        elseif (tmp2 < 1)
            TrajOut=cat(1,TrajOut,[ky,kz]); % The peripheral K-space sector(between inner and outer circle

        end

    end
end
for ky = double(ystart-5: yend+4)
    for kz = double(zstart-5:zend+4)
        %idx = (ky + ycenter) * KzN + (kz + zcenter);
        % 椭圆方程
        tmp1 = ky * ky / power(a*1.08,2) + kz * kz / power(b*1.08,2);
        % 椭圆面积S=π(圆周率)×a×b
        if (tmp1 < 1)% 得到椭圆内部从中心开始的25%的点
            TrajPeri=cat(1,TrajPeri,[ky,kz]);
        end
    end
end
InComp=complex(TrajIn(:,1),TrajIn(:,2));
OutComp=complex(TrajOut(:,1),TrajOut(:,2));
PeriComp=complex(TrajPeri(:,1),TrajPeri(:,2));
PeriComp=setdiff(PeriComp,InComp);
PeriComp=setdiff(PeriComp,OutComp);
TrajPeri=[real(PeriComp),imag(PeriComp)];

%% 测试各个区域有无重合点
TrajOutComplex=complex(TrajOut(:,1),TrajOut(:,2));
TrajInComplex=complex(TrajIn(:,1),TrajIn(:,2));
TrajPeriComplex=complex(TrajPeri(:,1),TrajPeri(:,2));

if ~isempty(intersect(TrajInComplex,TrajOutComplex))
    error('内外分配存在重合');
end
if ~isempty(intersect(TrajPeriComplex,TrajOutComplex))
    error('边缘与外环分配存在重合');
end

%% 调整内外点数

[TrajIn,TrajOut,TrajPeri]=pointAdjust(TrajIn,TrajOut,TrajPeri,NE,a2_2,b2_2,a1_2,b1_2);
TrajInComplex=complex(TrajIn(:,1),TrajIn(:,2));
TrajOutComplex=complex(TrajOut(:,1),TrajOut(:,2));
TrajPeriComplex=complex(TrajPeri(:,1),TrajPeri(:,2));

if ~isempty(intersect(TrajInComplex,TrajOutComplex))
    error('内外分配存在重合');
end
if ~isempty(intersect(TrajPeriComplex,TrajOutComplex))
    error('边缘与外环分配存在重合');
end
% figure,scatter(TrajIn(:,1),TrajIn(:,2),'b*')
% hold on
% scatter(TrajOut(:,1),TrajOut(:,2),'r+')
% hold on
% scatter(TrajPeri(:,1),TrajPeri(:,2),'b*')
nShotIn=size(TrajIn,1)/NE;
nShotOut=size(TrajOut,1)/NE;
nShotOutAve=nShotOut/4;
weight = abs(TrajIn(:, 1)) + 1e-6 * TrajIn(:, 1) + 1e-12 * TrajIn(:, 2);
% 使用 sort 进行排序
[~, idx] = sort(weight);
sortedTrajIn = TrajIn(idx, :);
%areaSan=floor(size(TrajOut,1)/NE);
%% 获取位置后分割外圈环型
%posSet=zeros(2,8,NE);
delta=[pi/2,pi/2,1.5*pi,1.5*pi];
LabelOut=[];
deltaSet=[];
%pointSet=[];
%TrajOutComplex=complex(TrajOut(:,1),TrajOut(:,2));

L=a;
ForeNumSet = zeros(4,1);
NextNumSet=ForeNumSet;
n_=hanning(NE);
%LastNuFore=0;
pointAllocated=[];
LastNumFore=0;
lambda=0;
T=zeros(4,1);
T_=zeros(4,1);
LoopCounters = cell(NE, 4,1);
LoopCounters{1,1}='k,k_,LastNum';
tic
for n=0:NE-2


    segTraj=sortedTrajIn(n*nShotIn+1:n*nShotIn+nShotIn,:);
    yPMax = max(segTraj(segTraj(:,1)>=0,1));%Ky正轴
    yPMin = min(segTraj(segTraj(:,1)>=0,1));

    yNMax = max(segTraj(segTraj(:,1)<0,1));%Ky负轴
    yNMin = min(segTraj(segTraj(:,1)<0,1));
%     if n==48
%         toc
%         tic
%     end

    %两个都不为空，即两边的都有数据
    % 获取正轴边缘点坐标
    [posPUpSet,posPBotSet]=getPeriPos(yPMin:yPMax,sortedTrajIn);
    [posNUpSet,posNBotSet]=getPeriPos(yNMin:yNMax,sortedTrajIn);
    posSetLine=[];
    posSetLine=cat(1,posSetLine,posPUpSet(end,:));%第一象限
    posSetLine=cat(1,posSetLine,posNUpSet(1,:));%第二象限
    posSetLine=cat(1,posSetLine,posNBotSet(1,:));%第三象限
    posSetLine=cat(1,posSetLine,posPBotSet(end,:));%第四象限

    theta=1/NE*1/3*pi/power(n+1,1/4);
    
    %% 获取外环扇形区域
    
    
    

    RepNum=0;
    %L=a*(1.5+3*n_(n+1));
    PosOutSet=[];
    NumSet=zeros(4,1);
    for xx=1:4
        if T(xx,1)==1
            NumSet(xx,1)=ForeNumSet(xx,1);
            NextNumSet(xx,1)=ForeNumSet(xx,1)-nShotOutAve;
            T_(xx,1)=1;
            continue
        else
            Num=0;
            k=0;
            k_=0;
            CountSet=[];
            while (Num+ForeNumSet(xx,1))<nShotOutAve
                k=k+1;
                [PosOut,ang,k_,lambda]=getArea(xx,n+1,NE,k_,lambda,nShotOutAve,Num,ForeNumSet(xx,1),TrajOut,L,posSetLine(xx,:),...
                    theta,delta(xx),pointAllocated,1);
                %[PosOut,RepNum]=RemoveDopPoint(PosOut,pointAllocated);
                Num=size(PosOut,1);
                Num_=Num;
                delta(xx)=ang;
                if Num+ForeNumSet(xx,1)>2*nShotOutAve
                    while Num+ForeNumSet(xx,1)>2*nShotOutAve
                        [PosOut,ang,k_,lambda]=getArea(xx,n+1,NE,k_,lambda,nShotOutAve,Num,ForeNumSet(xx,1),TrajOut,L,posSetLine(xx,:),...
                            theta,delta(xx),pointAllocated,0);
                        [PosOut,RepNum]=RemoveDopPoint(PosOut,pointAllocated);
                        Num=size(PosOut,1);
                        lam=lambda*power(RepNum,1/5)/5;
                        delta(xx)=ang;
                    end
                end
            %% 记录每次循环的次数
            Count=[k,k_,Num+ForeNumSet(xx,1)-nShotOutAve];
            CountSet=[CountSet;Count];
            
            end
            LoopCounters{n+2,xx}=CountSet;
            PosOut=RemoveDopPoint(PosOut,pointAllocated);
            PosOutSet=[PosOutSet;PosOut];
            NumSet(xx,1)=Num;
            NextNumSet(xx,1)=Num+ForeNumSet(xx,1)-nShotOutAve;
            T_(xx,1)=0;

        end
    end


    if n>0
        %% 检查两个相邻之间是否存在空缺点
        RegionIn=sortedTrajIn;%是否存在问题，应不应该算上上述的
        RegionOut=cat(1,TrajOut(LabelOut(:,1),:),[real(PosOutSet),imag(PosOutSet)]);
        RegionCondi=RegionOut(RegionOut(:,1)>=0 & RegionOut(:,2)>=0,:);
        Region=[RegionIn;RegionOut];
        [VaPointSet_,NumVaSet]=vacancyRemove(Region,RegionCondi,TrajPeri,TrajOutComplex);%分象限检查更快
        PosOutSet=[PosOutSet;VaPointSet_];
        PosOutSet=RemoveDopPoint(PosOutSet,pointAllocated);
        NumSet=NumSet+NumVaSet;
        NextNumSet=NextNumSet+NumVaSet;
        pointAllocated=unique(cat(1,pointAllocated,PosOutSet),'stable');
    else
        pointAllocated=cat(1,pointAllocated,PosOutSet);
    end


    %% 分配区域

    [PosForSelf,PosForNext,PosForNext2,T,~]=allocateOverPoint(PosOutSet,delta,posSetLine,NumSet+ForeNumSet-nShotOutAve,nShotOutAve,T_);
    if ~T_(1,1)
        LabelOut=AllocateLabel(PosForSelf.p1,PosForNext.p1,PosForNext2.p1,LabelOut,TrajOutComplex,n,T(1,1));
    end
    if ~T_(2,1)
        LabelOut=AllocateLabel(PosForSelf.p2,PosForNext.p2,PosForNext2.p2,LabelOut,TrajOutComplex,n,T(2,1));
    end
    if ~T_(3,1)
        LabelOut=AllocateLabel(PosForSelf.p3,PosForNext.p3,PosForNext2.p3,LabelOut,TrajOutComplex,n,T(3,1));
    end
    if ~T_(4,1)
        LabelOut=AllocateLabel(PosForSelf.p4,PosForNext.p4,PosForNext2.p4,LabelOut,TrajOutComplex,n,T(4,1));
    end
    if size(unique(LabelOut(:,1)))<size(LabelOut(:,1))
        error('存在重复序号');
    end
    %LastNumSet=size(LabelOut2,1);


    %% 初始值更新过为下一次循环做准备

    ForeNumSet=NextNumSet;
    %delta=[angSet(1,1),angSet(1,2),angSet(1,3),angSet(1,4)];
    deltaSet=cat(1,deltaSet,delta);
   
end
toc

numAll=1:size(TrajOutComplex,1);
idxLastAre=setdiff(numAll,LabelOut(:,1));
LabelOut=cat(1,LabelOut,[idxLastAre',repmat(NE,[size(idxLastAre,2) 1])]);



figure,scatter(sortedTrajIn(1:nShotIn,1),sortedTrajIn(1:nShotIn,2),'b*')
hold on
for i=1:floor((NE-1)/2)

    scatter(sortedTrajIn((2*i-1)*nShotIn+1:2*i*nShotIn,1),sortedTrajIn((2*i-1)*nShotIn+1:2*i*nShotIn,2),'r+')

    scatter(sortedTrajIn(2*i*nShotIn+1:(2*i+1)*nShotIn,1),sortedTrajIn(2*i*nShotIn+1:(2*i+1)*nShotIn,2),'b*')
end
scatter(sortedTrajIn((2*i+1)*nShotIn+1:NE*nShotIn,1),sortedTrajIn((2*i+1)*nShotIn+1:NE*nShotIn,2),'r+')

for i=1:NE/2

    scatter(TrajOut(LabelOut(LabelOut(:,2)==2*i-1),1),TrajOut(LabelOut(LabelOut(:,2)==2*i-1),2),'b*')
    scatter(TrajOut(LabelOut(LabelOut(:,2)==2*i),1),TrajOut(LabelOut(LabelOut(:,2)==2*i),2),'r+')
end

scatter(TrajOut(LabelOut(LabelOut(:,2)==NE),1),TrajOut(LabelOut(LabelOut(:,2)==NE),2),'g+')
title('Ky=180,Kz=40,NE=50,CenterRatio=25%');
if size(LabelOut,1) > size(TrajOut,1)
    warning('存在重复点');
end

sortedTrajOut=zeros(size(LabelOut,1),3);
numSet=zeros(NE,1);
for i=1:NE
    numSet(i)=size(find(LabelOut(:,2)==i),1);
    sortedTrajOut((i-1)*nShotOut+1:(i)*nShotOut,:)=[TrajOut(LabelOut(:,2)==i,:),repmat(i,[numSet(i) 1])];
    newSortedTrajIn((i-1)*nShotIn+1:(i)*nShotIn,:)=[sortedTrajIn((i-1)*nShotIn+1:(i)*nShotIn,:),repmat(i,[nShotIn 1])];
end
newTrajIn=ResortPoint(newSortedTrajIn,nShotIn,NE);
newTrajOut=ResortPoint(sortedTrajOut,nShotOut,NE);
% data=[newTrajIn;newTrajOut];
% filename = 'D:\\Matlab代码\\Traj_NE50.txt'; % 文件名






data=ones(1,KyN*KzN)*-10;


%% function section

function ResortedTraj=ResortPoint(Traj,n,NE)
ResortedTraj=Traj;
for i=1:n
    for j=1:NE
        ResortedTraj((i-1)*NE+j,:)=Traj((j-1)*n+i,:);
    end
end

end

function [TrajInNew,TrajOutNew,TrajPeri]= pointAdjust(TrajIn,TrajOut,TrajPeri,NE,a2_2,b2_2,a1_2,b1_2)
if mod(size(TrajIn,1),NE)~=0
    distances = power(TrajOut(:,1),2)/a1_2+power(TrajOut(:,2),2)/b1_2;

    % 根据距离和额外条件进行排序
    [~, indices] = sortrows(distances);
    % 根据排序结果重新排列点的顺序
    sortedTrajOut = TrajOut(indices,:);
    remainsIn=mod(size(TrajIn,1),NE);
    TrajInNew=cat(1,TrajIn,sortedTrajOut(1:NE-remainsIn,:));
    TrajOutNew=sortedTrajOut(NE-remainsIn+1:end,:);
else
    TrajInNew=TrajIn;
    TrajOutNew=TrajOut;
end
nShot=[size(findXxPoint(TrajOut,1),1),size(findXxPoint(TrajOut,2),1),size(findXxPoint(TrajOut,3),1),size(findXxPoint(TrajOut,4),1)]/NE;
nShotOutAve=round(mean(nShot));
[TrajOut1,TrajPeri1]=outAdjust(findXxPoint(TrajOutNew,1),findXxPoint(TrajPeri,1),NE,nShotOutAve,a2_2,b2_2);
[TrajOut2,TrajPeri2]=outAdjust(findXxPoint(TrajOutNew,2),findXxPoint(TrajPeri,2),NE,nShotOutAve,a2_2,b2_2);
[TrajOut3,TrajPeri3]=outAdjust(findXxPoint(TrajOutNew,3),findXxPoint(TrajPeri,3),NE,nShotOutAve,a2_2,b2_2);
[TrajOut4,TrajPeri4]=outAdjust(findXxPoint(TrajOutNew,4),findXxPoint(TrajPeri,4),NE,nShotOutAve,a2_2,b2_2);
TrajOutNew=[TrajOut1;TrajOut2;TrajOut3;TrajOut4];
TrajPeri=[TrajPeri1;TrajPeri2;TrajPeri3;TrajPeri4];
if mod(size(TrajInNew,1),NE)~=0 || mod(size(TrajOutNew,1),NE*4)~=0
    error('调整后的内圈或圆环点数不为回波数的整数');
end
end
function PointOut=findXxPoint(Point,xx)
switch xx
    case 1
        PointOut=Point(Point(:,1)>=0 & Point(:,2)>=0,:);

    case 2
        PointOut=Point(Point(:,1)<0 & Point(:,2)>=0,:);
    case 3
        PointOut=Point(Point(:,1)<0 & Point(:,2)<0,:);
    case 4
        PointOut=Point(Point(:,1)>=0 & Point(:,2)<0,:);
end
end
function [TrajOut,TrajPeri]=outAdjust(TrajOutIn,TrajPeriIn,NE,nShotOutAve,a2_2,b2_2)
if size(TrajOutIn,1)-nShotOutAve*NE~=0

    distances1 = power(TrajPeriIn(:,1),2)/a2_2+power(TrajPeriIn(:,2),2)/b2_2;
    distances2 = power(TrajOutIn(:,1),2)/a2_2+power(TrajOutIn(:,2),2)/b2_2;
    % 根据距离和额外条件进行排序
    [~, indices1] = sortrows(distances1);
    [~, indices2] = sortrows(distances2);
    % 根据排序结果重新排列点的顺序
    sortedTrajPeri = TrajPeriIn(indices1,:);
    sortedTrajOut = TrajOutIn(indices2,:);
    if size(TrajOutIn,1)<NE*nShotOutAve
        TrajOut=cat(1,sortedTrajOut,sortedTrajPeri(1:NE*nShotOutAve-size(TrajOutIn,1),:));
        TrajPeri=sortedTrajPeri(NE*nShotOutAve-size(TrajOutIn,1)+1:end,:);
    elseif size(TrajOutIn,1)>NE*nShotOutAve
        TrajOut=sortedTrajOut(1:NE*nShotOutAve,:);
        TrajPeri=cat(1,sortedTrajOut(NE*nShotOutAve+1:end,:),sortedTrajPeri);
    end

else
    TrajOut=TrajOutIn;
    TrajPeri=TrajPeriIn;
end


end
function [posUpSet,posBotSet]=getPeriPos(yRange,Traj)
posUpSet=[];
posBotSet=[];
for i=yRange
    zUpSet=max(Traj(Traj(:,1)==i,2));
    zBotSet=min(Traj(Traj(:,1)==i,2));

    posUpSet=cat(1,posUpSet,[i,zUpSet]);
    posBotSet=cat(1,posBotSet,[i,zBotSet]);
end
end

function LabelOut=AllocateLabel(PosForSelf,PosForNext,PosForNext2,LabelIn,TrajOutComplex,n,t)

[~, indices1] = ismember(PosForSelf, TrajOutComplex);
LabelOut1=cat(2,indices1,repmat(n+1,[size(indices1,1) 1]));
[~, indices2] = ismember(PosForNext, TrajOutComplex);
LabelOut2=cat(2,indices2,repmat(n+2,[size(indices2,1) 1]));
if size(LabelOut2,1)==0
    LabelOut2=[];
end

LabelOut=cat(1,LabelIn,LabelOut1,LabelOut2);
if t
    [~, indices3] = ismember(PosForNext2, TrajOutComplex);
    LabelOut3=cat(2,indices3,repmat(n+3,[size(indices3,1) 1]));
    if size(LabelOut3,1)==0
        LabelOut3=[];
    end
    LabelOut=cat(1,LabelOut,LabelOut3);
end

end
% 交点求解函数
function [PosForSelf,PosForNext,PosForNext2,T,NextNumSet]=allocateOverPoint(PosOutSet,delta,posSetLine,NumForNext,NshotAve,T_)
SortedPosOutPUp=sortPointsOnLine(PosOutSet(real(PosOutSet)>=0 & imag(PosOutSet)>=0),delta(1),posSetLine(1,:));
SortedPosOutNUp=sortPointsOnLine(PosOutSet(real(PosOutSet)<0 & imag(PosOutSet)>=0),delta(2),posSetLine(2,:));
SortedPosOutNBot=sortPointsOnLine(PosOutSet(real(PosOutSet)<0 & imag(PosOutSet)<0),delta(3),posSetLine(3,:));
SortedPosOutPBot=sortPointsOnLine(PosOutSet(real(PosOutSet)>=0 & imag(PosOutSet)<0),delta(4),posSetLine(4,:));
PosForSelf=struct();
PosForNext=struct();
PosForNext2=struct();
[PosForSelf.p1,PosForNext.p1,PosForNext2.p1,NextNum1]=allocateByCondi(SortedPosOutPUp,NumForNext(1,1),NshotAve,T_(1,1));
[PosForSelf.p2,PosForNext.p2,PosForNext2.p2,NextNum2]=allocateByCondi(SortedPosOutNUp,NumForNext(2,1),NshotAve,T_(2,1));
[PosForSelf.p3,PosForNext.p3,PosForNext2.p3,NextNum3]=allocateByCondi(SortedPosOutNBot,NumForNext(3,1),NshotAve,T_(3,1));
[PosForSelf.p4,PosForNext.p4,PosForNext2.p4,NextNum4]=allocateByCondi(SortedPosOutPBot,NumForNext(4,1),NshotAve,T_(4,1));
T=[size(PosForNext.p1,1)==NshotAve;size(PosForNext.p2,1)==NshotAve;size(PosForNext.p3,1)==NshotAve;size(PosForNext.p4,1)==NshotAve];
NextNumSet=[NextNum1;NextNum2;NextNum3;NextNum4];
end
function [PosForSelf,PosForNext,PosForNext2,NextNum]=allocateByCondi(SortedPosOutPUp,Num,NumAve,t)
if Num==0 || t==1
    PosForSelf=SortedPosOutPUp;PosForNext=[];PosForNext2=[];
    NextNum=0;
elseif Num<=NumAve
    PosForNext2=[];
    PosForNext=SortedPosOutPUp(1:Num);
    PosForSelf=SortedPosOutPUp(Num+1:end);
    NextNum=Num;
else
    PosForNext2=SortedPosOutPUp(1:Num-NumAve);
    PosForNext=SortedPosOutPUp(Num-NumAve+1:Num);
    PosForSelf=SortedPosOutPUp(Num+1:end);
    NextNum=Num-NumAve;
end
end
function sortedPoints = sortPointsOnLine(pointsComplex, ang_rad, Pos)
%points=cat(2,real(pointsComplex),imag(pointsComplex));
% 计算每个点到给定直线的距离
a=Pos(1);
b=Pos(2);
%
if b>0

    distances = abs(imag(pointsComplex) - (b-1) - tan(ang_rad) * (real(pointsComplex) - a))';

else
    distances = abs(imag(pointsComplex) - (b+1) - tan(ang_rad) * (real(pointsComplex) - a))';
end

abs_y = abs(imag(pointsComplex))';


extra_condition = (-abs_y);


[~, indices] = sortrows([distances; extra_condition]', [1 -2]);

sortedPoints = pointsComplex(indices);

end


function PosVa=Dect4conn(PosAreaPart)
PosAreaComp=complex(PosAreaPart(:,1),PosAreaPart(:,2));
xSet=sort(unique(PosAreaPart(:,1)));
PosVa=[];
for i=1:size(xSet,1)-2
    x=xSet(i+1);
    yMax=max(PosAreaPart(PosAreaPart(:,1)==x,2));
    yMin=min(PosAreaPart(PosAreaPart(:,1)==x,2));
    condition=ones(1,5);
    for y=yMin+1:yMax-1
        condition(1) = ~ismember(complex(x,y), PosAreaComp);
        if ~condition(1)
            continue;
        else
            condition(2) = ismember(complex(x,y+1), PosAreaComp);
            condition(3) = ismember(complex(x,y-1), PosAreaComp);
            condition(4) = ismember(complex(x-1,y), PosAreaComp);
            condition(5)= ismember(complex(x+1,y), PosAreaComp);
        end
        if(all(condition) )
            PosVa=cat(1,PosVa,complex(x,y));
        end
    end
end
end

function [VaPointSet,NumSet]=vacancyRemove(RegionInOut,RegionCondi,TrajPeri,TrajOutComplex)
%第一象限
Region1=RegionInOut(RegionInOut(:,1)>=0 & RegionInOut(:,2)>=0,:);
Peri1=TrajPeri(TrajPeri(:,1)>=0 & TrajPeri(:,2)>=0,:);

%第二象限
Region2=RegionInOut(RegionInOut(:,1)<0 & RegionInOut(:,2)>=0,:);
Peri2=TrajPeri(TrajPeri(:,1)<0 & TrajPeri(:,2)>=0,:);


%第三象限
Region3=RegionInOut(RegionInOut(:,1)<0 & RegionInOut(:,2)<0,:);
Peri3=TrajPeri(TrajPeri(:,1)<0 & TrajPeri(:,2)<0,:);

%第四象限
Region4=RegionInOut(RegionInOut(:,1)>=0 & RegionInOut(:,2)<0,:);
Peri4=TrajPeri(TrajPeri(:,1)>=0 & TrajPeri(:,2)<0,:);


pointsVa1=Dect4conn([Region1;Peri1]);
pointsVa2=Dect4conn([Region2;Peri2]);
pointsVa3=Dect4conn([Region3;Peri3]);
pointsVa4=Dect4conn([Region4;Peri4]);

pointsVa1=intersect(pointsVa1,TrajOutComplex);
pointsVa2=intersect(pointsVa2,TrajOutComplex);
pointsVa3=intersect(pointsVa3,TrajOutComplex);
pointsVa4=intersect(pointsVa4,TrajOutComplex);
%VaPointSet=struct();

%边界检查

if  min(RegionCondi(:,2))==0
    %    if loge==1
    pointsBound=Dect4conn([RegionInOut;TrajPeri]);

    %     else
    %         pointsBound=Dect4conn(RegionInOut);
    %
    %     end
    pointsVa1=unique([pointsVa1;pointsBound(real(pointsBound)>=0 & imag(pointsBound)>=0)]);
    pointsVa2=unique([pointsVa2;pointsBound(real(pointsBound)<0 & imag(pointsBound)>=0)]);
    pointsVa3=unique([pointsVa3;pointsBound(real(pointsBound)<0 & imag(pointsBound)<0)]);
    pointsVa4=unique([pointsVa4;pointsBound(real(pointsBound)>=0 & imag(pointsBound)<0)]);
end
VaPointSet=[pointsVa1;pointsVa2;pointsVa3;pointsVa4];
NumSet=[size(pointsVa1,1);size(pointsVa2,1);size(pointsVa3,1);size(pointsVa4,1)];

end

function [PosOutUpdata,RepNum]=RemoveDopPoint(PosOut,pointAllocated)


PosOutUpdata=setdiff(PosOut,pointAllocated);

RepNum=size(PosOut,1)-size(PosOutUpdata,1);
end
function [PosOut,ang_,w,lambda]=getArea(xx,n,NE,w,lambda,nShotOut,Num,LastNum,TrajOut,L,pos,...
    theta,ang,pointAllocated,loge)
if loge==1
    w=0;
    while(nShotOut > Num+LastNum)
        w=w+1;
        dif=(nShotOut-Num-LastNum);
        %conditions=[angPUp(2)<0-0.01,angNUp(1)>pi+0.01,angNBot(1)<pi-0.01,angPBot(2)>2*pi+0.01];
        lambda=2*power(dif,1/2)/power(n,1/4);
        
        switch xx
            case {1,3}
                ang=ang-theta*lambda;
                ang_=ang;
            case {2,4}
                ang=ang+theta*lambda;
                ang_=ang;
        end
        PosSetor = findSectorPoints(pos, ang, L,n,NE, xx);
        PosOut = unique(intersect(complex(TrajOut(:,1),TrajOut(:,2)),PosSetor));
        [PosOut,~]=RemoveDopPoint(PosOut,pointAllocated);
        Num=size(PosOut,1);


    end
elseif loge==0
    w=w+1;
    dif=(Num+LastNum-nShotOut);
    lambda_=lambda/(power(w,1/4))/(power(n,1/4));
    switch xx
        case {1,3}
            ang=ang+theta*lambda_;
            ang_=ang;
        case {2,4}
            ang=ang-theta*lambda_;
            ang_=ang;
    end
    PosSetor = findSectorPoints(pos, ang, L,n,NE, xx);
    PosOut = unique(intersect(complex(TrajOut(:,1),TrajOut(:,2)),PosSetor));
    [PosOut,~]=RemoveDopPoint(PosOut,pointAllocated);
    Num_=size(PosOut,1);

else
    error('增减面积的标志输入有误');

end
%PosPeri = intersect(complex(TrajPeri(:,1),TrajPeri(:,2)),PosSetor);
end
function PosFinal = findSectorPoints(posLine, ang, L,n,NE, LabelQuad)
% 找到满足条件的扇形的整型坐标点
% 参数:
% a, b: 第一条直线过的点的坐标
% c, d: 第二条直线过的点的坐标
% k1: 第一条直线的斜率
% k2: 第二条直线的斜率
% L: 扇形的边界长度
% direction: 扇形的方向，可以选择 1 或 -1，1 表示 a,b 到交点的方向，-1 表示交点到 a,b 的方向
a=posLine(1,1);
b=posLine(1,2);
[x0, y0] = findIntersection(a, b, ang);


%lambda=n-1;
gamma=pi/(4+(NE-20)/10);
ang_=gamma;
% 根据方向确定起始角度和结束角度
switch LabelQuad
    case 1
        angle_start = ang+ang_;
        angle_end = ang;
    case 2
        angle_start = ang-ang_;
        angle_end = ang;
    case 3
        angle_start = ang+ang_;
        angle_end = ang;
    case 4
        angle_start = ang-ang_;
        angle_end = ang;
end
% 生成扇形的边界点坐标
angNum=round(20*power(n,1/4));
%angNum=400;
%angShift=(angle_end-angle_start)/angNum*2;
theta_ = linspace(angle_start, angle_end, angNum);
LSet=linspace(1, L,L)';%改成500可能提高运行速度
% 计算扇形内部的坐标点
x_float = x0 + LSet*cos(theta_);
y_float = y0 + LSet*sin(theta_);

% 转换成整型坐标

x_list=round(x_float(:));
y_list=round(y_float(:));
complex_points = complex(x_list, y_list);
%unique_Pos = unique(complex_points);
unique_Pos = complex_points;
PosAll=cat(2,real(unique_Pos),imag(unique_Pos));
Pos=findLabelQuadrantPoints(PosAll,LabelQuad);
PosFinal = complex(Pos(:,1), Pos(:,2));
end
function [x, y] = findIntersection(a, b,  theta)

    % 计算直线的斜率
    slope = tan(theta);
    
    % 计算与y轴的交点
    x = 0;
    y = b - a * slope;
end
function Pos = findLabelQuadrantPoints(PosAll,LabelQuad)
switch LabelQuad
    case 1
        logIdx=(PosAll(:,1)>=0) & (PosAll(:,2)>=0);
        Pos=PosAll(logIdx,:);

    case 2
        logIdx=(PosAll(:,1)<0) & (PosAll(:,2)>=0);
        Pos=PosAll(logIdx,:);

    case 3
        logIdx=(PosAll(:,1)<0) & (PosAll(:,2)<0);
        Pos=PosAll(logIdx,:);

    case 4
        logIdx=(PosAll(:,1)>=0) & (PosAll(:,2)<0);
        Pos=PosAll(logIdx,:);

    otherwise
        error('象限标志输入有误，必须输入1、2、3或4代表象限的标志符')
end

end



