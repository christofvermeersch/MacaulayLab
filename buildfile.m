function plan = buildfile(~)
    %BUILDFILE   Tasks to run before building the project.
    %   plan = BUILDFILE(~) creates a plan for MATLAB's buildtool.

    % Copyright (c) 2024 - Christof Vermeersch

    % Import built-in tasks:
    import matlab.buildtool.tasks.CodeIssuesTask
    import matlab.buildtool.tasks.TestTask

    % Create a plan from local task functions:
    plan = buildplan(localfunctions);

    % Add two default tasks (code issues and tests) to the plan:
    plan("check") = CodeIssuesTask;
    plan("test") = TestTask;
    
    % Make the 'archive' task dependent on the two default tasks:
    plan("archive").Dependencies = ["check" "test"];

    % Make the 'archive' task the default task in the plan:
    plan.DefaultTasks = "archive";
    
end

function archiveTask(~)
    %ARCHIVETASK   Creation of an archive.
    %   ARCHIVETASK(~) creates a .zip archive of the toolbox.
    filename = "macaulaylab_" + ...
        string(datetime("today",Format="yyyymmdd"));
    zip(filename,"*")
end