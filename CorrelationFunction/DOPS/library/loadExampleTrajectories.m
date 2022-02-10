function loadExampleTrajectories(NameOfFile,Path2Project,reducedDataSet)
    
    % set default value if necessary
    if nargin > 2
        loadReducedDataSet=reducedDataSet;
    else
        loadReducedDataSet=false;
    end
    
    % set file path for dataset
    if loadReducedDataSet
         Path2File=[ Path2Project NameOfFile '_reduced.mat' ];
    else
         Path2File=[ Path2Project NameOfFile '.mat' ];
    end
     
    % load data
    Data=load(Path2File); 

    if loadReducedDataSet
        Data=Data(1).DataRed;
    end
     
    Lipid_H=Data.(NameOfFile);
    ps=1e-12;
    %-------------------------------
    DeltaT=1.0*ps;     % time steps
    %-------------------------------

    Traj_x=squeeze(Lipid_H(:,1,:));
    Traj_y=squeeze(Lipid_H(:,2,:));
    Traj_z=squeeze(Lipid_H(:,3,:));


    %% -------------Clean workspace
    clearvars -except   DeltaT ps Traj_x Traj_y Traj_z Mol Path2Save
    
end