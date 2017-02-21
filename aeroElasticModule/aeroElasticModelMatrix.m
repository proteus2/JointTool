classdef aeroElasticModelMatrix < aeroElasticModel
    %LINSTATICAEROSTRUCTURALMODEL Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (SetAccess = protected)
        Tas
        Tsa
        pos
    end
    properties (Access = private)
        Ras
        das
    end
    
    methods
        function obj = aeroElasticModelMatrix(inAeroModel,inStrModel) %constructor
           obj = obj@aeroElasticModel(inAeroModel,inStrModel); %call constructor of superclass
           % Coupling matrices
           % A --> S (Nearest Neighbour): on force application point
           [obj.Tas,obj.Ras] = obj.a2s(obj.aeroModel.fvap',obj.strModel.grid);
           % S --> A (Nearest Neighbour): on aerodynamic grid points
           [obj.Tsa,~,obj.pos] = obj.a2s(obj.aeroModel.grid',obj.strModel.grid);
           % Arms for moment calculation
           obj.das = obj.Ras;
        end
        function obj = transformForces(obj)
            Mx =  obj.aeroModel.F(2,:).*obj.das(:,3)' - obj.aeroModel.F(3,:).*obj.das(:,2)';
            My =  obj.aeroModel.F(1,:).*obj.das(:,3)' - obj.aeroModel.F(3,:).*obj.das(:,1)';
            Mz =  obj.aeroModel.F(1,:).*obj.das(:,2)' - obj.aeroModel.F(2,:).*obj.das(:,1)';
            Fs = obj.Tas*[obj.aeroModel.F; Mx; My; Mz]';
%             keyboard
            obj.strModel.Fs = reshape(Fs',size(Fs,1)*size(Fs,2),1);
        end
        function obj = transformDisplacements(obj)
            % Recover rotations and upgrade grid
            obj.aeroModel.gridDeflected = obj.rotateGrid()';
        end
    end
    methods
        function [Tas,Ras,posv] = a2s(obj,aeronodes, strnodes)
            VisualCheck = 0;

            na = size(aeronodes,1);
            ns = size(strnodes,1);

            Tas  = zeros(ns,na);
            Ras  = zeros(ns,3);
            posv = zeros(ns,1); 

            for i = 1:na
               Xas        = [strnodes(:,1) - aeronodes(i,1), strnodes(:,2) - aeronodes(i,2), strnodes(:,3) - aeronodes(i,3)];
               % Absolute Distance
               R          = (Xas(:,1).^2 + Xas(:,2).^2 + Xas(:,3).^2).^(1/2);
               [~,idx]    = min(R);

                % I (no relaxation scheme)
            %    Tas(pos,i) = 1;
            %    Ras(i,:)   = -Xas(pos,:);

                % II (3 nodes relaxation scheme)
               if idx~=1 && idx~=ns
                  Tas(idx-1,i) = Tas(idx-1,i) + 1/3;
                  Tas(idx  ,i) = Tas(idx  ,i) + 1/3;
                  Tas(idx+1,i) = Tas(idx+1,i) + 1/3;
                  Ras(i,:)     = -1/3*(Xas(idx-1,:) + Xas(idx,:) + Xas(idx+1,:));
               elseif idx==1
                  Tas(idx  ,i) = Tas(idx  ,i) + 1/2; 
                  Tas(idx+1,i) = Tas(idx+1,i) + 1/2; 
                  Ras(i,:)     = -1/2*(Xas(idx,:) + Xas(idx+1,:));
               elseif idx==ns
                  Tas(idx-1,i) = Tas(idx-1,i) + 1/2; 
                  Tas(idx  ,i) = Tas(idx  ,i) + 1/2; 
                  Ras(i,:)     = -1/2*(Xas(idx-1,:) + Xas(idx,:));
               end

               % Store node indices
               posv(i,1) = idx;
            end

            if VisualCheck
                figure()
                hold on
                plot3(aeronodes(:,1),aeronodes(:,2),aeronodes(:,3),'k.')
                plot3(strnodes(:,1),strnodes(:,2),strnodes(:,3),'ro-')
                colorid = repmat([1 ; 2],length(strnodes),1);
                col{1} = 'go'; col{2} = 'yo';
                for i=1:length(strnodes(:,1))
            %         figure()
            %         hold on
            %         plot3(aeronodes(:,1),aeronodes(:,2),aeronodes(:,3),'k.')
            %         plot3(strnodes(:,1),strnodes(:,2),strnodes(:,3),'ro-')
                    for j=1:length(Tas(i,:))
                        if Tas(i,j)~=0
                            plot3(aeronodes(j,1),aeronodes(j,2),aeronodes(j,3),col{colorid(i)})
            %                 plot3(aeronodes(j,1),aeronodes(j,2),aeronodes(j,3),'go')
                        end
                    end
                end
                axis equal
            end
        end
        function varargout = rotateGrid(obj)
            VisualCheck = 0;

            disp = reshape(obj.strModel.disp,6,size(obj.strModel.disp,1)/6)';
            T = disp(:,1:3);
            R = disp(:,4:6);

            Agrid = obj.aeroModel.geometry.grid';
            Sgrid = obj.strModel.grid;

            Tsa = obj.Tsa';

            mapT = Tsa*T;
            mapR = Tsa*R;

            pos = obj.pos;

            Parot = zeros(size(Tsa,1),3);
            for i=1:size(Tsa,1)
                RotMat = expon(mapR(i,:));
                Pa = Agrid(i,:);
                Pel = Sgrid(pos(i),:);
                Parot(i,:) = RotMat*(Pa - Pel)' + Pel' + mapT(i,:)';
            end

            varargout{1} = Parot;

            if VisualCheck
                figure()
                hold all
                plot3(Agrid(:,1),Agrid(:,2),Agrid(:,3),'rx')
                plot3(Parot(:,1),Parot(:,2),Parot(:,3),'ko')
                axis equal
            end

            function [R,dR,d2R]=expon(t,dt)

            if nargin==1
                flag = 1;
                dR = 0;
            else
                flag = 2;
            end

            I  = eye(3,3);

            al  = sqrt(t(1)^2+t(2)^2+t(3)^2);

            Rsk = skew(t);

            if al<eps
                R  = I;
            else
                Rsk2        = Rsk^2;
                R  = I+sin(al)/al*Rsk+0.5*(sin(al/2)/(al/2))^2*Rsk2;
            end

            if flag==2
                % dRdx
                [Rsk,dRsk]  = skew(t,dt);
                if al<eps
                    dR = dRsk;
                else
                    dal = 1/2/sqrt(t(1)^2+t(2)^2+t(3)^2)*(2*t(1).*dt(1,:)+2*t(2).*dt(2,:)+2*t(3).*dt(3,:));
                    c1 = 1/al^2*(cos(al)*al-sin(al))*dal;
                    c2 = 4*sin(al/2)/al^3*(.5*cos(al/2)*al-sin(al/2))*dal;
                    for i=1:size(Rsk,1)
                        for j=1:size(Rsk,2)
                            Rskc1(i,j,:) = c1.*Rsk(i,j);
                            Rskc2(i,j,:) = c2.*Rsk2(i,j);
                        end    
                    end
                    for i=1:size(dt,2)
                        RdRdRR(:,:,i) = Rsk*dRsk(:,:,i)+dRsk(:,:,i)*Rsk;
                    end
                    dR = Rskc1+sin(al)/al.*dRsk+Rskc2+2*(sin(al/2)/al)^2.*RdRdRR;
                end

                %d2Rdx2
                if al<eps
                    d2R = zeros(3,3,size(dt,2),size(dt,2));
                else
                    for i=1:size(dt,2)
                        for j=1:size(dt,2)
                            d2al = -1/al^3*(t(1).*dt(1,j)+t(2).*dt(2,j)+t(3).*dt(3,j))*(t(1).*dt(1,i)+t(2).*dt(2,i)+t(3).*dt(3,i))+...
                                1/al*(dt(1,j)*dt(1,i)+dt(2,j)*dt(2,i)+dt(3,j)*dt(3,i));

                            c1 = 1/al^2*(cos(al)*al-sin(al));
                            c2 = (al^2*(-al*sin(al))-(cos(al)*al-sin(al))*2*al)/al^4;
                            c3 = 4*(al^3*1/2*cos(al/2)-sin(al/2)*3*al^2)/al^6*(al/2*cos(al/2)-sin(al/2));
                            c4 = 4*sin(al/2)/al^3*(-1/2*al/2*sin(al/2));
                            c5 = 4*sin(al/2)/al^3*(.5*cos(al/2)*al-sin(al/2));
                            c6 = 0.5*(sin(al/2)/(al/2))^2;

                            d2R(:,:,i,j) = c1*Rsk*d2al+...
                                c1*dRsk(:,:,j)*dal(:,i)+...
                                c2*Rsk*dal(:,i)*dal(:,j)+...
                                c1*dal(:,j)*dRsk(:,:,i)+...
                                c3*Rsk2*dal(:,i)*dal(:,j)+...
                                c4*Rsk2*dal(:,i)*dal(:,j)+...
                                c5*dal(:,i)*(Rsk*dRsk(:,:,j)+dRsk(:,:,j)*Rsk)+...
                                c5*Rsk2*d2al+...
                                c5*dal(:,j)*(Rsk*dRsk(:,:,i)+dRsk(:,:,i)*Rsk)+...
                                c6*(dRsk(:,:,i)*dRsk(:,:,j)+dRsk(:,:,j)*dRsk(:,:,i));
                        end
                    end
                end
            end
            end 
        end
    end
end
