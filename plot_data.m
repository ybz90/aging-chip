function plot_data

    % READ FROM EXT COFNIG FILE S.T. FOR EVERY CELL WE KNOW WHAT ITS LIFE
    % IS FRAME RANGE AS MANUALLY COUNTED
    % WHEN THAT CELL IS BEING PLOTTED, IT'S DATA WILL BE 1:LIFE instead of
    % 1:sz(1), which is the full frame #

  % Yuan Zhao 05/07/2015

  positions = 27:29; %read this from configfile


  % Array for storing all trajectory data across all cells
  all_traj = [];

  % Horizontally concatenate traj matrices for every position, forming a super array with dimensions:
  % # frames x # cells/traps (from all positions) x # fluorescent channel
  for i = positions
    pos = num2str(i);
    traj_file = ['xy',pos,'/xy',pos,'_traj.mat']
    load(traj_file);
    all_traj = horzcat(all_traj,traj);
    sz = size(all_traj); % get the dimensions of the combined all_traj super array
  end


  % INIT SUPER PLOTS FOR EACH FLU


  % For every cell in all_traj...
  for j = 1:sz(2)

      colors = ['r','g','b']; %styles eventually not just colors?
      %figure;
      intensities = [];
      %X = 1:sz(1);
      X = 1:400;
      % SET TO EXACT LIFE BASED ON MANUAL CHECKING; WE CAN MAP TO THAT $ BY
      % LOOKING AT POSITION XY AND 1-7 CELL/COLN; ORIGI ID IS EZ, JUST
      % TAKE J IN THE CONCAT MAT AND DIV BY 7 AKA COLN; REMAINDER IS THE
      % CELL # 1-7 AND QUOTIENT INTEGER IS THE POSITION # PAST FIRST XY

      for flu = 1:sz(3) %for each flu channel

          intensity = all_traj(X,j,flu);

          intensities{flu} = intensity;

      end

      B = cell2mat(intensities(1));
      C = cell2mat(intensities(2));

      [ax,p1,p2] = plotyy(X,B,X,C,'plot');
      %grid(ax(1),'on')
      p1.LineWidth = 2;
      p2.LineWidth = 2;
      p2.LineSTyle = '--';
      ylabel(ax(1),'GFP');
      ylabel(ax(2),'iRFP nuc');
      xlabel(ax(2),'time in frames');
      title('POSITION CELL # AND ALSO LIFESPAN');

      % ADD TO MEGA PLOT FOR EACH FLU
  end
end


  % PLOT ALL TOGETHER ON ONE GRAPH IN ADDIDITION TO EACH CELL SEPARATELY
