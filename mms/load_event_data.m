%event = '20150815_1';


switch event    
    case '20150815_1'
        cd /Users/Cecilia/Data/MMS/2015Aug15/
        tsl = 7364480;

        FID = fopen('mms3_scb_20150815_125000t.dat');
        Bscmt = fread(FID,tsl,'int64');
        Bscmt = int64(Bscmt);
        Bscmt = EpochTT(Bscmt);

        FID = fopen('mms3_scb_20150815_125000Bxyz.dat');
        B = fread(FID,tsl*3,'float');
        B = [B(1:tsl) B((1+tsl):2*tsl) B((1+2*tsl):3*tsl)];

        B3sc = TSeries(Bscmt,B,'to',1);

        %load brst mode interval
        if 0 % load sburt electric field and survey magnetic field data
        tmpDataObj = dataobj('data/mms3_edp_brst_ql_dce2d_20150815130004_v0.2.0.cdf');
        Exyz = get_variable(tmpDataObj,'mms3_edp_dce_xyz_dsl');
        Exyz = mms.variable2ts(Exyz);
        Exyz.data(find(abs(Exyz.data) > 200)) = NaN;

        tmpDataObj = dataobj('data/mms3_dfg_srvy_ql_20150815_v0.0.3.cdf');
        Bxyz = get_variable(tmpDataObj,'mms3_dfg_srvy_dmpa');
        Bxyz = mms.variable2ts(Bxyz);

        tlimit = irf.tint(Exyz.time.start.utc,Exyz.time.stop.utc);
        Exyz = Exyz.tlim(tlimit);
        Bscm = Bscm.tlim(tlimit);
        tlimitl = tlimit+[-60 60];
        Bxyz = Bxyz.tlim(tlimitl);
        Bxyz = Bxyz.resample(Exyz);
        Bscm = Bscm.resample(Exyz);
        end
    case '20150828_1'
        cd /Users/Cecilia/Data/MMS/2015Aug28/
        if 1 % load C3 
            % Search coil data
            tsl = 401280;            
            FID = fopen('mms3_scb_20150828_150600t.dat');
            Bscmt = fread(FID,tsl,'int64');
            Bscmt = int64(Bscmt);
            Bscmt =  EpochTT(Bscmt);
            FID = fopen('mms3_scb_20150828_150600Bxyz.dat');
            Bscm = fread(FID,tsl*3,'float');
            Bscm = [Bscm(1:tsl) Bscm((1+tsl):2*tsl) Bscm((1+2*tsl):3*tsl)];
            B3sc = TSeries(Bscmt,Bscm,'to',1);

            % Electric field data
            tmpDataObj = dataobj('data/mms3_edp_brst_ql_dce2d_20150828150624_v0.2.0.cdf');
            dslE3 = get_variable(tmpDataObj,'mms3_edp_dce_xyz_dsl');
            dslE3 = mms.variable2ts(dslE3);
            dslE3.data(find(abs(dslE3.data) > 200)) = NaN;
            
            % Survey magnetic field data
            tmpDataObj = dataobj('data/mms3_dfg_srvy_ql_20150828_v0.0.3.cdf');
            B3fg = get_variable(tmpDataObj,'mms3_dfg_srvy_dmpa');
            B3fg = mms.variable2ts(B3fg);

            % Limit burst data to search coil interval
            tlimit = irf.tint(B3sc.time.start.utc,B3sc.time.stop.utc);
            dslE3 = dslE3.tlim(tlimit);
            B3sc = B3sc.tlim(tlimit);
            tlimitl = tlimit+[-60 60];
            B3fg = B3fg.tlim(tlimitl);
            B3fg = B3fg.resample(dslE3);
        end
        if 1 % load C4
            % Search coil data
            tsl = 729600;
            FID = fopen('mms4_scb_20150828_125300t.dat');
            Bscmt = fread(FID,tsl,'int64');
            Bscmt = int64(Bscmt);
            Bscmt =  EpochTT(Bscmt);

            FID = fopen('mms4_scb_20150828_125300Bxyz.dat');
            Bscm = fread(FID,tsl*3,'float');
            Bscm = [Bscm(1:tsl) Bscm((1+tsl):2*tsl) Bscm((1+2*tsl):3*tsl)];
            B4sc = TSeries(Bscmt,Bscm,'to',1);

            % Electric field burst data
            tmpDataObj = dataobj('data/mms4_edp_brst_ql_dce2d_20150828125314_v0.2.0.cdf');
            dslE4 = get_variable(tmpDataObj,'mms4_edp_dce_xyz_dsl');
            dslE4 = mms.variable2ts(dslE4);
            dslE4.data(find(abs(dslE4.data) > 200)) = NaN;

            tmpDataObj = dataobj('data/mms4_dfg_srvy_ql_20150828_v0.0.3.cdf');
            B4fg = get_variable(tmpDataObj,'mms4_dfg_srvy_dmpa');
            B4fg = mms.variable2ts(B4fg);
            
            % Limit burst data to search coil interval      
            tlimit = irf.tint(B4sc.time.start.utc,B4sc.time.stop.utc);
            dslE4 = dslE4.tlim(tlimit);
            B4sc = B4sc.tlim(tlimit);
            tlimitl = tlimit+[-60 60];
            B4fg = B4fg.tlim(tlimitl);
            B4fg = B4fg.resample(dslE4);
        end
        
        
end