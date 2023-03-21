close all, clear all
disp('Running...')

%% Configure optical path
tilt_ang = 0.05;
f = [50,100,50,75,150];
x = [0,100,250,320,445,670,820];                                            % [start, lens_pos_1, lens_pos 2, ..., focal plane]
lens_dia = 100;
light_type = 'parallel';                                                    % parallel and diverse

%% Initial optical path
nSteps = length(f) + length(x);
x_ = zeros(1,nSteps);
d_ = zeros(1,nSteps);
f_ = zeros(1,nSteps);
x_idx = 1;
f_idx = 1;
for nStep = 1:nSteps
    if nStep ==1
        x_(nStep) = x(1);
        d_(nStep) = 0;
        f_(nStep) = 0;
        x_idx = x_idx+1;
    elseif nStep == nSteps
        x_(nStep) = x(end);
        d_(nStep) = x(end)-x(end-1);
        f_(nStep) = 0;
    elseif mod(nStep,2)==0
        x_(nStep) = x(x_idx);
        d_(nStep) = x(x_idx)-x(x_idx-1);
        f_(nStep) = 0;
    else
        x_(nStep) = x(x_idx);
        d_(nStep) = 0;
        f_(nStep) = f(f_idx);
        x_idx = x_idx+1;
        f_idx = f_idx+1;
    end
end

%% Simulate optical path
figure, set(gcf, 'windowstate','maximized'), title('Perfect alignment')
for y0 = -10:2.5:10
    switch light_type
        case 'parallel'
            u = [y0;tilt_ang];
        case 'diverse'
            if y0>0
                u = [y0;tilt_ang];
            elseif y0==0
                u = [y0;0];
            else
                u = [y0;-tilt_ang];
            end
    end
    
    % calculate optical matrix
    y(1) = u(1);
    for nStep = 2:length(x_)
        % transport operator
        if (d_(nStep)~=0)
            P = [1, d_(nStep); 0, 1];
            u = P*u;
        end
        
        % lens operator
        if (f_(nStep)~=0)
            hold on, line([x_(nStep),x_(nStep)],[-lens_dia/2,lens_dia/2],'color','k','linewidth',2)
            text(x_(nStep)+5,5,['f = ',num2str(f_(nStep))],'BackgroundColor','w')
            L = [1, 0; -1/f_(nStep), 1];
            u = L*u;
        end
        
        y(nStep) = u(1);
    end
    hold on, plot(x_,y,'r')
end
hold on, line([x_(end),x_(end)],[-100,100],'color','c','linewidth',2)
text(x_(nStep)+5,5,'focal plane','BackgroundColor','w')
xticks([min(x):25:max(xlim)])
daspect([1,1,1]), grid on

%%
disp('Done.')

