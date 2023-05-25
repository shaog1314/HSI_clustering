function [flag, updateChange, numEqualInARowNew, stoppingCause] = stoppingCriterionUV(iii,U,U2,V,V2,param,numEqualInARow)
%Evaluates if the stopping criterion of the used method is fulfilled or
%not based on the relative change of the matrices.

%INPUT:
%iii:               Name of Dataset
%U:                 New cluster membership matrix iteration
%U2:                Old cluster membership matrix iteration
%V:                 New centroid matrix iteration
%V2:                Old centroid matrix iteration
%param:             Struct, which contains all predefined hyperparameters.
%numEqualInARow:    Number of times that the cluster membershipt
%                   assignments has not changed yet.

%OUTPUT:
%flag:              Stopping flag
%updateChange:      Struct; Relative change of all matrix iterates
%numequalInARowNew: The number of times, how often there were no changes in
%                   a row in the cluster membership assignments
%stoppingCause:     Cause of setting the flag to true (string)

% Written by Pascal Fernsel
% (Center for Industrial Mathematics, University of Bremen,
% p.fernsel@uni-bremen.de)

% Reference paper: 
% P. Fernsel, "Spatially Coherent Clustering Based on Orthogonal
% Nonnegative Matrix Factorization", Journal of Imaging, 2021.

% This code comes with no guarantee or warranty of any kind.

    updateChange.U = norm(U-U2, 'fro')/norm(U2, 'fro');
    updateChange.V = norm(V-V2, 'fro')/norm(V2, 'fro');
    
    stoppingCause='';
    
    switch param.stoppingCrit
        case 'relChange'
            flag = false;
            if iii >= param.maxIt
                flag = true;
                stoppingCause = 'maxIt';
            elseif mod(iii,2) == 0
                if updateChange.U < param.stopLimUV && ...
                        updateChange.V < param.stopLimUV
                  flag = true;
                  stoppingCause = 'relChange';
                end
            end
            numEqualInARowNew = numEqualInARow;
        case 'clusterAssign'
            flag = false;
            UHardNew=convertToHardClustering(U);
            UHardOld=convertToHardClustering(U2);


            %Permute here first of all the columns, in case there are just some
            %changes in the order of the features.
            for i=1:size(UHardNew,2)
                sameColIndex=find(all(~(UHardNew-UHardOld(:,i)),1));
                if length(sameColIndex)==1
            %             colPermuted=true;
                    swapIndex=1:size(UHardNew,2);
                    swapIndex(i)=sameColIndex;
                    swapIndex(sameColIndex)=i;
                    UHardNew=UHardNew(:, swapIndex);
                elseif length(sameColIndex)>1
                    error('There are equal columns in the hard clustering of U.')
                end
            end

           %Check for equality (except permutation of columns)
            if all(all(UHardNew==UHardOld))
               numEqualInARowNew = numEqualInARow + 1;
               if numEqualInARowNew>=param.numEqualInARowMax
                   flag=true;
                   stoppingCause = 'clusterAssign';
               end
            else
               numEqualInARowNew=0;
            end

            if iii >= param.maxIt
                flag = true;
                stoppingCause = 'maxIt';
            end
        case 'maxIt'
            flag = false;
            if iii >= param.maxIt
                flag = true;
                stoppingCause = 'maxIt';
            end
            numEqualInARowNew = numEqualInARow;
        otherwise
            error('Unknown choice of the stopping criterion.\n')
    end
    
    
    
end