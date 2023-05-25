function [kappa, AA,OA,ACC,NMI,ARI,PUR,grp] = fast_performance_overview(gnd,idx)
%CALC_ACC Function that calculates the accuracy of a clustering algorithm.
% Inputs: gnd - Nx1 vector of ground truth labels
%         idx - Nx1 vector of estimated labels
% Output: kappa coefficient, AA: averaged accuracy, OA: overall accuracy,
% ACC: accuracy per class, grp: final clustering results
n1=find(gnd~=0);
gnd_1=gnd(n1);
idx_1=idx(n1);
best_c = bestMap_1(gnd_1,idx_1);


Label1 = unique(gnd_1);
Label2 = unique(idx_1);
nClass2 = length(Label2);
grp = zeros(size(idx));
for i=1:nClass2
    grp(idx == Label2(i)) = Label1(best_c(i));
end


grp1=grp;
grp1(gnd==0)=0;
[kappa, AA,OA,ACC, PUR] = evaluate_results(grp1(:), gnd(:));
NMI = getNMI(gnd(n1), grp(n1));
ARI = rand_index(gnd(n1), grp(n1), 'adjusted');


    function [kappa, acc_A,acc_O,acc,pur] = evaluate_results(gtest, grd)
        if sum(grd==0)~=0
            c=max(grd)-min(grd);
            CM1=confusionmat(gtest,double(grd)); % g_test  vs grd 
            CM=CM1(2:end,2:end);
            d=sum(diag(CM));
            dCM=diag(CM);%pixels are classified brilliantly
            Ci=sum(CM1,1);
            Cj=sum(CM1,2);
            Call=sum(Ci)-CM1(1,1);
            acc_O=d/Call;
            kappa=(Call*d-Ci(2:end)*Cj(2:end))/(Call^2-Ci(2:end)*Cj(2:end));
            
            %groud truth
            acc=zeros(c,1);
            for i1=1:c
                acc(i1)=dCM(i1)/Ci(i1+1);
            end
            acc_A=sum(acc)/c;
            pur=sum(max(CM,[],2))/Call;
        else
            c=max(grd)-min(grd)+1;
            CM=confusionmat(gtest,double(grd));
            d=sum(diag(CM));
            dCM=diag(CM);
            Ci=sum(CM,1);
            Cj=sum(CM,2);
            Call=sum(Ci);
            acc_O=d/Call;
            kappa=(Call*d-Ci*Cj)/(Call^2-Ci*Cj);
            
            %groud truth
            acc=zeros(c,1);
            for i1=1:c
                acc(i1)=dCM(i1)/Ci(i1);
            end
            acc_A=sum(acc)/c;
            pur=sum(max(CM,[],2))/Call;
        end 
    end
end