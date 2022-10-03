classdef PolarCode_SCL < handle
    % Polar code
    properties
        
        block_length;
        info_length;
        n;
        info_set;
        A;
        
        %scl decoder
        list_size;
        p_scl;
        c_scl;
        llr_path_metric;
        
        inactivePathIndices;
        inactivePathIndicesSize;
        activePathArray;
        
        % CRC
        crc_matrix;
        crc_size;
        
    end
    
    methods
        %% constructor
        function obj = PolarCode_SCL(block_length, info_length, info_set, crc_size, crc_matrix)
            
            obj.block_length = block_length;
            obj.info_length = info_length;
            obj.n = log2(block_length);
            obj.info_set = info_set;  % the number of ones in info_set is expected to be equal to (info_length + crc_size)
            obj.A = find(obj.info_set);
            
            if nargin < 4
                crc_size = 0;
                crc_matrix = randi([0 1], obj.crc_size, info_length);
            end
            obj.crc_size = crc_size;
            obj.crc_matrix = crc_matrix;
                        
        end
        
        
        %%  Decoder helper
        
        function [b] = crc_check(obj, info_bits)
            crc_check =  mod(obj.crc_matrix*info_bits, 2)';
            b = all(crc_check == 0);
        end
        
        
        %%  SCL decoder
        
        function [u, x, llrout] = decode_scl_llr(obj, llrin, list_size)
            
            obj.list_size = list_size;
            obj.initializeDataStructures();
            l_index = obj.assignInitialPath();
            obj.p_scl{l_index + 1}(:,1) = llrin;
            obj.polar_decode_scl(obj.info_set, 0, 0);
            
            l_index_max = obj.findMostProbablePath(1);
            u = obj.c_scl{l_index_max + 1}(:, obj.n + 1);
            x = obj.c_scl{l_index_max + 1}(:, 1);
            llrout = obj.p_scl{l_index_max + 1}(:, obj.n + 1);
            
        end
        
        function [u, x, llrout, path_metric_list] = decode_scl_llr_getAllList(obj, llrin, list_size)
            
            obj.list_size = list_size;
            obj.initializeDataStructures();
            l_index = obj.assignInitialPath();
            obj.p_scl{l_index + 1}(:,1) = llrin;
            obj.polar_decode_scl(obj.info_set, 0, 0);
            
            u = cell(1, obj.list_size);
            x = cell(1, obj.list_size);
            llrout = cell(1, obj.list_size);
            for l_index = 0 : obj.list_size - 1
                u{l_index + 1} = obj.c_scl{l_index + 1}(:, obj.n + 1);
                x{l_index + 1} = obj.c_scl{l_index + 1}(:, 1);
                llrout{l_index + 1} = obj.p_scl{l_index + 1}(:, obj.n + 1);
            end
            path_metric_list = obj.llr_path_metric;
            
        end
        
        function polar_decode_scl(obj, info_set, phi, lambda)
            
            N = length(info_set);
            rate0_pattern = zeros(1,N);
            if (N==1)
                if info_set == 0
                    obj.continuePaths_FrozenBit(phi);
                else
                    obj.continuePaths_UnfrozenBit(phi);
                end
            elseif all(info_set == rate0_pattern)
                obj.continuePaths_Rate0Node(N, lambda, phi);
            else
                info_set_1 = info_set(1:N/2);
                info_set_2 = info_set(N/2+1:end);
                for l_index = 0 : obj.list_size - 1
                    if obj.activePathArray(l_index + 1) == 0
                        continue;
                    end
                    obj.p_scl{l_index + 1}(phi+1:phi+N/2, lambda+2) = PolarCode_SCL.cnop_llr_minsum(obj.p_scl{l_index + 1}(phi+1:phi+N/2, lambda+1), obj.p_scl{l_index + 1}(phi+N/2+1:phi+N, lambda+1));
                end
                obj.polar_decode_scl(info_set_1, phi, lambda + 1);
                for l_index = 0 : obj.list_size - 1
                    if obj.activePathArray(l_index + 1) == 0
                        continue;
                    end
                    obj.p_scl{l_index + 1}(phi+N/2+1:phi+N, lambda+2) = PolarCode_SCL.vnop_llr(obj.p_scl{l_index + 1}(phi+1:phi+N/2, lambda+1), obj.p_scl{l_index + 1}(phi+N/2+1:phi+N, lambda+1), ...
                        obj.c_scl{l_index + 1}(phi+1:phi+N/2, lambda+2));
                end
                obj.polar_decode_scl(info_set_2, phi+N/2, lambda + 1);
                for l_index = 0 : obj.list_size - 1
                    if obj.activePathArray(l_index + 1) == 0
                        continue;
                    end
                    obj.c_scl{l_index + 1}(phi+1:phi+N/2, lambda+1) = mod(obj.c_scl{l_index + 1}(phi+1:phi+N/2, lambda+2) + obj.c_scl{l_index + 1}(phi+N/2+1:phi+N, lambda+2), 2);
                    obj.c_scl{l_index + 1}(phi+N/2+1:phi+N, lambda+1) = obj.c_scl{l_index + 1}(phi+N/2+1:phi+N, lambda+2);
                end
            end
            
        end
                
        function initializeDataStructures(obj)
            
            obj.inactivePathIndices = zeros(obj.list_size,1);
            obj.inactivePathIndicesSize = 0;
            % the above two variables are used to define a stack
            
            obj.activePathArray =  zeros(obj.list_size,1);
            
            obj.llr_path_metric =  zeros(obj.list_size, 1);
            obj.p_scl = cell(1, obj.list_size);
            obj.c_scl = cell(1, obj.list_size);
            for i_list = 1:obj.list_size
                obj.p_scl{i_list} = zeros(obj.block_length, obj.n + 1);
                obj.c_scl{i_list} = zeros(obj.block_length, obj.n + 1);
            end
            
            for i_list = 0 : obj.list_size - 1
                obj.activePathArray(i_list + 1) = 0;
                obj.inactivePathIndices(i_list + 1) = i_list;
            end
            
            obj.inactivePathIndicesSize  = obj.list_size;
            
        end
        
        function l_index = assignInitialPath(obj)
            l_index = obj.inactivePathIndices(obj.inactivePathIndicesSize);
            obj.inactivePathIndicesSize = obj.inactivePathIndicesSize - 1;
            obj.activePathArray(l_index + 1) = 1;
        end
        
        function l_index_clone = clonePath(obj, l_index)
            l_index_clone = obj.inactivePathIndices(obj.inactivePathIndicesSize);
            obj.inactivePathIndicesSize = obj.inactivePathIndicesSize - 1;
            obj.activePathArray(l_index_clone + 1) = 1;
            
            obj.llr_path_metric(l_index_clone + 1) = obj.llr_path_metric(l_index + 1);
            obj.p_scl{l_index_clone + 1} = obj.p_scl{l_index + 1};
            obj.c_scl{l_index_clone + 1} = obj.c_scl{l_index + 1};
        end
        
        function killPath(obj, l_index)
            obj.activePathArray(l_index + 1) = 0;
            obj.inactivePathIndices(obj.inactivePathIndicesSize + 1) = l_index;
            obj.inactivePathIndicesSize = obj.inactivePathIndicesSize + 1;
            obj.llr_path_metric(l_index + 1) = 0;
        end
        
        function continuePaths_FrozenBit(obj, phi)
            
            for l_index = 0 : obj.list_size -1
                
                if obj.activePathArray(l_index + 1) == 0
                    continue;
                end
                
                obj.c_scl{l_index + 1}(phi + 1, obj.n + 1) = 0;
                
                obj.llr_path_metric(l_index + 1) = obj.llr_path_metric(l_index + 1) + log1p(exp(-obj.p_scl{l_index + 1}(phi + 1, obj.n + 1)));
                
            end
            
        end
                
        function continuePaths_UnfrozenBit(obj, phi)
            
            probForks = realmax * ones(obj.list_size, 2);
            index = 0;
            for l_index = 0 : obj.list_size - 1
                
                if obj.activePathArray(l_index + 1)
                    probForks(l_index + 1, 1) =  obj.llr_path_metric(l_index + 1) + log1p(exp(-obj.p_scl{l_index + 1}(phi + 1, obj.n + 1)));
                    probForks(l_index + 1, 2) =  obj.llr_path_metric(l_index + 1) + log1p(exp(obj.p_scl{l_index + 1}(phi + 1, obj.n + 1)));
                    index = index + 1;
                end
                
            end
            
            rho = min(2*index, obj.list_size);
            contForks = zeros(obj.list_size, 2);
            prob = sort(probForks(:));
            threshold = prob(rho);
            
            num_populated = 0;
            for l_index = obj.list_size-1 : -1 : 0
                for j_index = 1 : 2
                    if num_populated == rho
                        break;
                    end
                    if  probForks(l_index + 1, j_index) <= threshold
                        contForks(l_index + 1, j_index) = 1;
                        num_populated = num_populated + 1;
                    end
                end
            end
            
            
            for l_index = 0 : obj.list_size - 1
                if obj.activePathArray(l_index + 1) == 0
                    continue;
                end
                
                if (contForks(l_index + 1, 1) == 0) && (contForks(l_index + 1, 2) == 0)
                    obj.killPath(l_index);
                end
            end
            
            for l_index = 0 : obj.list_size - 1
                if obj.activePathArray(l_index + 1) == 0
                    continue;
                end
                if (contForks(l_index + 1, 1) == 0) && (contForks(l_index + 1, 2) == 0)
                    continue;
                end
                if (contForks(l_index + 1, 1) == 1) && (contForks(l_index + 1, 2) == 1)
                    l_index_clone = obj.clonePath(l_index);
                    
                    if obj.p_scl{l_index + 1}(phi + 1, obj.n + 1) > 0
                        obj.c_scl{l_index + 1}(phi + 1, obj.n + 1) = 0;
                        obj.c_scl{l_index_clone + 1}(phi + 1, obj.n + 1) = 1;

                        obj.llr_path_metric(l_index + 1) = obj.llr_path_metric(l_index + 1) + log1p(exp(-obj.p_scl{l_index + 1}(phi + 1, obj.n + 1)));
                        obj.llr_path_metric(l_index_clone + 1) = obj.llr_path_metric(l_index_clone + 1) + log1p(exp(obj.p_scl{l_index_clone + 1}(phi + 1, obj.n + 1)));
                    else
                        obj.c_scl{l_index + 1}(phi + 1, obj.n + 1) = 1;
                        obj.c_scl{l_index_clone + 1}(phi + 1, obj.n + 1) = 0;

                        obj.llr_path_metric(l_index + 1) = obj.llr_path_metric(l_index + 1) + log1p(exp(obj.p_scl{l_index + 1}(phi + 1, obj.n + 1)));
                        obj.llr_path_metric(l_index_clone + 1) = obj.llr_path_metric(l_index_clone + 1) + log1p(exp(-obj.p_scl{l_index_clone + 1}(phi + 1, obj.n + 1)));
                    end
                else
                    if contForks(l_index + 1, 1) == 1
                        obj.c_scl{l_index + 1}(phi + 1, obj.n + 1) = 0;
                        obj.llr_path_metric(l_index + 1) = obj.llr_path_metric(l_index + 1) + log1p(exp(-obj.p_scl{l_index + 1}(phi + 1, obj.n + 1)));
                    else
                        obj.c_scl{l_index + 1}(phi + 1, obj.n + 1) = 1;
                        obj.llr_path_metric(l_index + 1) = obj.llr_path_metric(l_index + 1) + log1p(exp(obj.p_scl{l_index + 1}(phi + 1, obj.n + 1)));
                    end
                end
                
            end
            
        end
        
        function continuePaths_Rate0Node(obj, node_size, lambda, phi)
            
            for l_index = 0 : obj.list_size -1
                
                if obj.activePathArray(l_index + 1) == 0
                    continue;
                end
                
                obj.llr_path_metric(l_index + 1) = obj.llr_path_metric(l_index + 1) + sum(log1p(exp(-obj.p_scl{l_index + 1}(phi + 1: phi + node_size, lambda + 1))));
                obj.c_scl{l_index + 1}(phi + 1: phi + node_size, lambda + 1:obj.n + 1) = 0;
                obj.p_scl{l_index + 1}(phi + 1: phi + node_size, lambda + 2:obj.n + 1) = obj.llr_path_metric(l_index + 1) * ones(node_size, obj.n - lambda);
                
            end
            
        end
        
        function l_index_max = findMostProbablePath(obj, crc_check)
            l_index_max = obj.list_size - 1;
            path_metric_min = realmax;
            
            path_with_crc = 0;
            for l_index = obj.list_size-1 : -1 : 0
                
                if obj.activePathArray(l_index + 1) == 0
                    continue;
                end
                
                if (crc_check) && (obj.crc_size ~= 0)
                    info_bits = obj.c_scl{l_index + 1}(obj.A, obj.n + 1);
                    if obj.crc_check(info_bits) == 0
                        continue;
                    end
                end
                
                path_with_crc = 1;
                if path_metric_min > obj.llr_path_metric(l_index + 1)
                    path_metric_min = obj.llr_path_metric(l_index + 1);
                    l_index_max = l_index;
                end
                
            end
            
            if (crc_check) && (path_with_crc == 0) % no path with crc check found
                l_index_max = obj.findMostProbablePath(0);
            end
            
        end
                
        
    end
    
    methods (Static)
        
        function z = cnop_llr(x, y)
            z = 2 * atanh(tanh(x/2) .* tanh(y/2));
        end
        
        function z = cnop_llr_minsum(x, y)
            xabs = abs(x);
            yabs = abs(y);
            maxi = max(xabs,yabs);
            mini = min(xabs,yabs);
            z = sign(x).*sign(y).*(mini+log1p(exp(-(xabs+yabs)))-log1p(exp(-(maxi-mini))));
        end
        
        function z = vnop_llr(x, y, u1)
            z = (1-2*u1) .* x + y;
        end
        
    end
    
end

