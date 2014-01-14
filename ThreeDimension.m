function ThreeDimension(action)
% ThreeDimension shoe how to control a 3D Ball with words.

global oDiv oDivMin colVec radius mcList d word wordNum ...
            aA lastA lastB aSize tspeed howElliptical dtr distr...
            beforeMouseX beforeMouseY currMouseX currMouseY

if nargin == 0
    action = 'start';
end

switch(action)
    
    %% initail & open Handle Graphics
    case 'start'
        
        oDiv = [300, 300];          % 200 x 200
        oDivMin = min(oDiv);
        radius = oDivMin * 0.8 / 2;
        dtr = pi/180;
        d = 300;

        colVec = rand(1,3);
        word = {'apple','boy','cat','dog','egg','fog','girl','heel',...
            'ip','jeep','kate','letter','mouse','nut','orange','pie'};
        wordNum = length(word);

        aA = zeros(wordNum,5);
        aA( : , 1) = 1 : wordNum;
        % offsetWidth, offsetHeight, cx, cy, cz, x, y, scale, alpha
        mcList = zeros( wordNum, 9);
        % active = 0;  % a flag
        lastA = 1;
        lastB = 1;
        distr = 1;
        tspeed = 10;
        aSize = 250;

        howElliptical = 1;

        for ind = 1:size(aA,1)
            mcList(ind,1:2) = [4 5];  % 
        end

        % 	initPositionAll;
        [mcList, aA] = initPositionAll( mcList, distr, aA);
        
        axis equal
        axis([0 oDiv(1) 0 oDiv(2)]);
        box on;     % add a frame
        title('3D Ball Show with Matlab');
        set(gcf, 'WindowButtonDownFcn', 'ThreeDimension down');
        
%         circle(radius, oDiv(1)./2, oDiv(2)./2);
        
        for ii = 1:wordNum
            hold on
            h = text(aA(ii,2), aA(ii,3), word(ii));
            fontsize = ceil(9 + 5*aA(ii,4)/radius);
            if (aA(ii,4)>=0)
                set(h, 'color', 'r','fontsize', fontsize);
            else
                set(h, 'color', 'black','fontsize', fontsize);
            end
        end
        hold off
        
        
        
    %% Callbacks when mouse down
    case 'down'
        clc
        % set Callbacks「tmouse move」 when mouse moves
        set(gcf, 'WindowButtonMotionFcn', 'ThreeDimension move');
        % set Callbacks「tmouse up」when mouse up
        set(gcf, 'WindowButtonUpFcn', 'ThreeDimension up');
        % print「Mouse down!」讯息
        fprintf('Mouse down!\n');
        
        currPt = get(gca, 'CurrentPoint');
        beforeMouseX = currPt(1,1);
        beforeMouseY = currPt(1,2);
        
        
    %% Callbacks when mouse moves
    case 'move'
        % 获得鼠标的当前位置
        currPt = get(gca, 'CurrentPoint');
        currMouseX = currPt(1,1);
        currMouseY = currPt(1,2);
%         centerX = oDiv(1) / 2;
%         centerY = oDiv(2) / 2;
        
        plot([currMouseX], [currMouseY],'*','color',colVec);
%         plot([centerX, currMouseX], [centerY, currMouseY],'-','color',colVec);
        colVec = mod(ceil(colVec*255 + 1),256)./256;
        
        axis equal
        axis([0 oDiv(1) 0 oDiv(2)]);
%         circle(radius, oDiv(1)./2, oDiv(2)./2);
        
        % print「Mouse is moving!」and its position
        fprintf('Mouse is moving! Current location = (%g, %g)\n', currMouseX, currMouseY);
        
        % 更新位置信息
        [mcList, lastA, lastB] = updatePos(mcList, currMouseX, beforeMouseX, currMouseY, beforeMouseY, lastA, lastB, aSize, radius, tspeed, howElliptical, d);
        aA = doPosition(mcList, oDiv, aA);
        
        for ii = 1:wordNum
            hold on
            h = text(aA(ii,2), aA(ii,3), word(ii));
            fontsize = ceil(9 + 5*aA(ii,4)/radius);
            if (aA(ii,4)>=0)
                set(h, 'color', 'r','fontsize', fontsize);
            else
                set(h, 'color', 'black','fontsize', fontsize);
            end
        end
        hold off
        
        beforeMouseX = currMouseX;
        beforeMouseY = currMouseY;
        
        %%   depthSort()
        
    %% Callbacks when mouse up
    case 'up'
        
        % clean move Callbacks
        set(gcf, 'WindowButtonMotionFcn', '');
        % clean up Callbacks
        set(gcf, 'WindowButtonUpFcn', '');
        % print「Mouse up!」
        fprintf('Mouse up!\n');
        
        
end

%% 
function [mcList,aA] = initPositionAll( mcList, distr, aA)

	max = size(mcList, 1);
	% 随机排列位置
    for i = 1 : max
        if ( distr == 1 )
			phi = acos( -1 + ( 2*i - 1 ) / max );
			theta = sqrt( max * pi ) * phi;
        else
			phi = rand * pi;
			theta = rand * 2 * pi;
        end
    % 坐标变换
        % cx
		mcList(i, 3) = radius * cos(theta) * sin(phi);
		% cy
        mcList(i, 4) = radius * sin(theta) * sin(phi);
		% cz
        mcList(i, 5) = radius * cos(phi);
		% left
		aA(i,2) = mcList(i, 3) + oDiv(1)/2 - mcList(i,1)/2 ;
		% top
        aA(i,3) = mcList(i, 4) + oDiv(2)/2 - mcList(i,2)/2;
        % ????
        aA(i,4) = mcList(i, 5);
    end
end 

%%
function aA = doPosition(mcList, oDiv,aA)
        l = oDiv(1,1) / 2;
        t = oDiv(1,2) / 2;
        for index = 1:size(mcList,1)
            % left
            aA(index, 2) = mcList(index,3) + l - mcList(index,1) / 2;
            % top
            aA(index, 3) = mcList(index,4) + t - mcList(index,2) / 2;
            %
            aA(index, 4) = mcList(index, 5);            
            % fontSize
            aA(index, 5) = ceil(12 * mcList(index,8) / 2)+8;
%             aA[i].style.filter="alpha(opacity="+100*mcList[i].alpha+")";
%             aA[i].style.opacity=mcList(i,9);
        end
end

%%
function [mcList, lastA, lastB] = updatePos(mcList, currMouseX, beforeMouseX, currMouseY, beforeMouseY, lastA, lastB, aSize, radius, tspeed, howElliptical, d)

%     if( 5 * mouseX > oDiv(1) || 5 * mouseY > oDiv(2) || mouseX < 0 || mouseY < 0)
%         acceleratorA = lastA * 0.98;
%         acceleratorB = lastB * 0.98;
%     else
%         acceleratorA = (min( max( -mouseY, -aSize ), aSize ) / radius ) * tspeed;
%         acceleratorB = (min( max( -mouseX, -aSize ), aSize ) / radius ) * tspeed;
%     end

    centerX = oDiv(1)/2;
    centerY = oDiv(2)/2;
    
    if (currMouseX-centerX)^2 + (currMouseY-centerY)^2 > radius^2
        speedX = 0;
        speedY = 0;
    else
        speedX = (currMouseX-beforeMouseX)/2;
        speedY = (currMouseY-beforeMouseY)/2;
        if(currMouseX-beforeMouseX<0 && currMouseY-beforeMouseY<0)
            speedX = -(currMouseX-beforeMouseX)/2;
            speedY = (currMouseY-beforeMouseY)/2;
        elseif (currMouseX-beforeMouseX>0 && currMouseY-beforeMouseY>0)
            speedX = -(currMouseX-beforeMouseX)/2;
            speedY = (currMouseY-beforeMouseY)/2;
        elseif (currMouseX-beforeMouseX>0 && currMouseY-beforeMouseY<0)
            speedX = (currMouseX-beforeMouseX)/2;
            speedY = -(currMouseY-beforeMouseY)/2;
        elseif (currMouseX-beforeMouseX<0 && currMouseY-beforeMouseY>0);
            speedX = (currMouseX-beforeMouseX)/2;
            speedY = -(currMouseY-beforeMouseY)/2;
        end
    end
    acceleratorA = speedX;
    acceleratorB = speedY;
    
    lastA = acceleratorA;
    lastB = acceleratorB;

    if(abs(acceleratorA) <= 0.01 && abs(acceleratorB) <= 0.01)
        return;
    end

    acceleratorC = 0;
    
    [ sa, ca, sb, cb, sc, cc ] = sineCosine(acceleratorA, acceleratorB, acceleratorC, dtr);
    
    % 
    for j = 1:wordNum
        
        rx1 = mcList(j,3);
        ry1 = mcList(j,4) * ca + mcList(j,5) * (-sa);
        rz1 = mcList(j,4) * sa + mcList(j,5) * ca;

        rx2 = rx1 * cb + rz1 * sb;
        ry2 = ry1;
        rz2 = rx1 * (-sb) + rz1 * cb;

        rx3 = rx2 * cc + ry2 * (-sc);
        ry3 = rx2 * sc + ry2 * cc;
        rz3 = rz2;

        mcList(j,3) = rx3;
        mcList(j,4) = ry3;
        mcList(j,5) = rz3;

        per = d / (d + rz3);

        mcList(j,6) = (howElliptical * rx3 * per) - ( howElliptical*2 );
        mcList(j,7) = ry3 * per;
        mcList(j,8) = per;
        mcList(j,9) = per;

        mcList(j,9) = ( mcList(j,9) - 0.6) * (10/6);

    end
end

%%
function [ sa, ca, sb, cb, sc, cc ] = sineCosine( a, b, c, dtr )

	sa = sin(a * dtr);
	ca = cos(a * dtr);
	sb = sin(b * dtr);
	cb = cos(b * dtr);
	sc = sin(c * dtr);
	cc = cos(c * dtr);
end


end
