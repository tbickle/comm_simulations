function [s_est] = mpam_map(r,M)
    % MAP detection estimation

    temp = zeros(log2(M),length(r));
    switch M
    case 2
%         temp = r-1;
        temp = (r+1)/2;
    case 4
%         for a = 1:length(r)
%             if (r(a) < -2)
%                 temp(:,a) = [0;0];
%             elseif (-2 <= r(a) && r(a) < 0)
%                 temp(:,a) = [1;0];
%             elseif (0 <= r(a) && r(a) < 2)
%                 temp(:,a) = [1;1];
%             elseif(2 <= r(a))
%                 temp(:,a) = [0;1];
%             end
%         end

%         for a = 1:length(r)
%             if (r(a) == 1)
%                 temp(:,a) = [0;0];
%             elseif (r(a) == 2)
%                 temp(:,a) = [1;0];
%             elseif (r(a) == 3)
%                 temp(:,a) = [1;1];
%             elseif(r(a) == 4)
%                 temp(:,a) = [0;1];
%             end
%         end

%         for a = 1:length(r)
%             if (r(a) == -3)
%                 temp(:,a) = [0;0];
%             elseif (r(a) == -1)
%                 temp(:,a) = [1;0];
%             elseif (r(a) == 1)
%                 temp(:,a) = [1;1];
%             elseif(r(a) == 3)
%                 temp(:,a) = [0;1];
%             end
%         end

        for a = 1:length(r)
            if (r(a) == -13)
                temp(:,a) = [0;0];
            elseif (r(a) == -4)
                temp(:,a) = [1;0];
            elseif (r(a) == 4)
                temp(:,a) = [0;1];
            elseif(r(a) == 13)
                temp(:,a) = [1;1];
            end
        end
    case 8
%         r = sign(r);
%         for a = 1:length(r)
%             f = r(1,a); g = r(2,a);
%             if (f==1 && g==0)
%                 temp(:,a) = [0;0;0];
%             elseif (f==1 && g==1)
%                 temp(:,a) = [1;0;0];
%             elseif (f==0 && g==1)
%                 temp(:,a) = [1;0;1];
%             elseif(f==-1 && g==1)
%                 temp(:,a) = [0;0;1];
%             elseif(f==-1 && g==0)
%                 temp(:,a) = [0;1;1];
%             elseif(f==-1 && g==-1)
%                 temp(:,a) = [1;1;1];
%             elseif(f==0 && g==-1)
%                 temp(:,a) = [1;1;0];
%             elseif(f==1 && g==-1)
%                 temp(:,a) = [0;1;0];
%             end
%         end
    case 16
%         r = round(r*10);
%         for a = 1:length(r)
%             f = r(1,a); g = r(2,a);
%             if (f==10 && g==0)       % 1
%                 temp(:,a) = [0;0;0;0];
%             elseif (f==9 && g==4)   % 2
%                 temp(:,a) = [1;0;0;0];
%             elseif (f==7 && g==7)   % 3
%                 temp(:,a) = [1;1;0;0];
%             elseif(f==4 && g==9)    % 4
%                 temp(:,a) = [0;1;0;0];
%             elseif(f==0 && g==10)    % 5
%                 temp(:,a) = [0;1;1;0];
%             elseif(f==-4 && g==9)   % 6
%                 temp(:,a) = [1;1;1;0];
%             elseif(f==-7 && g==7)   % 7
%                 temp(:,a) = [1;0;1;0];
%             elseif(f==-9 && g==4)   % 8
%                 temp(:,a) = [0;0;1;0];
%             elseif (f==-10 && g==0)  % 9
%                 temp(:,a) = [0;0;1;1];
%             elseif (f==-9 && g==-4) % 10
%                 temp(:,a) = [1;0;1;1];
%             elseif (f==-7 && g==-7) % 11
%                 temp(:,a) = [1;1;1;1];
%             elseif(f==-4 && g==-9)  % 12
%                 temp(:,a) = [0;1;1;1];
%             elseif(f==0 && g==-10)   % 13
%                 temp(:,a) = [0;1;0;1];
%             elseif(f==4 && g==-9)	% 14
%                 temp(:,a) = [1;1;0;1];
%             elseif(f==7 && g==-7)   % 15
%                 temp(:,a) = [1;0;0;1];
%             elseif(f==9 && g==-4)   % 16
%                 temp(:,a) = [0;0;0;1];
%             else                        % default
%                 temp(:,a) = [0;0;0;0];
%             end
%         end
    end

    s_est = temp;

end

