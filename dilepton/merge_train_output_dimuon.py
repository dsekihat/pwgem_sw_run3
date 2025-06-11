import subprocess

#__________________________________________________________
def merge(filename, original_tasknames, outfilename, outdirname):
    print("merge tasks in", filename);
    subprocess.run(["rootls", filename]);

    command_hadd = [];
    command_hadd.append("hadd");
    command_hadd.append("-f");
    command_hadd.append(outfilename);
    outnames_tmp = [];

    for original_taskname in original_tasknames:
        outfilename_tmp = filename.replace(".root", "_{0}.root".format(original_taskname));
        command_hadd.append(outfilename_tmp);
        outnames_tmp.append(outfilename_tmp);
        command = ["rootcp", "-r", "--recreate", "{0}:{1}".format(filename, original_taskname), "{0}:{1}".format(outfilename_tmp, outdirname)];
        print(command);
        subprocess.run(command);

    print(command_hadd);
    subprocess.run(command_hadd);

    for outname in outnames_tmp:
        command_rm = ["rm", "-rf", outname];
        print(command_rm);
        subprocess.run(command_rm);

#__________________________________________________________
def merge_occupancy(filename, taskname_prefix):
    list_centmin = [ 0, 10, 30, 50, 70,  0, 50];
    list_centmax = [10, 30, 50, 70, 90, 90, 90];
    list_occmin_src = [   0, 1000, 5000]; #source
    list_occmax_src = [1000, 5000, 99999]; #source

    list_occmin_dst = [   0]; #destination
    list_occmax_dst = [5000]; #destination
    list_taskname_suffix = [""]; #source

    ncen = len(list_centmin);
    nocc_src = len(list_occmin_src);
    nocc_dst = len(list_occmin_dst);
    for icen in range(0, ncen):
        cen1 = list_centmin[icen];
        cen2 = list_centmax[icen];
        for taskname_suffix in list_taskname_suffix:
            for iocc_dst in range(0, nocc_dst):
                occ1_dst = list_occmin_dst[iocc_dst];
                occ2_dst = list_occmax_dst[iocc_dst];
                print("desired occupancy {0:d} - {1:d}".format(occ1_dst, occ2_dst));
                list_original_tasknames = [];
                outfilename = filename.replace(".root", "_{0:02d}{1:02d}_occupancy{2:d}_{3:d}{4}.root".format(cen1, cen2, occ1_dst, occ2_dst, taskname_suffix));
                outdirname = "{0}_{1:02d}{2:02d}_occupancy{3:d}_{4:d}{5}".format(taskname_prefix, cen1, cen2, occ1_dst, occ2_dst, taskname_suffix);
                for iocc_src in range(0, nocc_src):
                    occ1_src = list_occmin_src[iocc_src];
                    occ2_src = list_occmax_src[iocc_src];

                    if (occ1_dst <= occ1_src and occ1_src <= occ2_dst) and (occ1_dst <= occ2_src and occ2_src <= occ2_dst):
                        taskname = "{0}_{1:02d}{2:02d}_occupancy{3:d}_{4:d}{5}".format(taskname_prefix, cen1, cen2, occ1_src, occ2_src, taskname_suffix);
                        list_original_tasknames.append(taskname);
                print(filename, list_original_tasknames, outfilename, outdirname);
                merge(filename, list_original_tasknames, outfilename, outdirname);

#__________________________________________________________
def merge_centrality(filename, taskname_prefix):
    list_centmin_src = [ 0, 10, 30, 50, 70]; #source
    list_centmax_src = [10, 30, 50, 70, 90]; #source

    list_centmin_dst = [ 0]; #destination
    list_centmax_dst = [90]; #destination

    list_occmin = [""]; #source
    list_occmax = [""]; #source

    list_taskname_suffix = ["", "_newsel8"]; #source

    ncen_src = len(list_centmin_src);
    ncen_dst = len(list_centmin_dst);
    nocc = len(list_occmin);
    for iocc in range(0, nocc):
        occ1 = list_occmin[iocc];
        occ2 = list_occmax[iocc];
        for taskname_suffix in list_taskname_suffix:
            for icen_dst in range(0, ncen_dst):
                cen1_dst = list_centmin_dst[icen_dst];
                cen2_dst = list_centmax_dst[icen_dst];
                print("desired centrality {0:d} - {1:d}".format(cen1_dst, cen2_dst));
                list_original_tasknames = [];

                if occ1 == "" and occ2 == "":
                    outfilename = filename.replace(".root", "_{0:02d}{1:02d}{2}.root".format(cen1_dst, cen2_dst, taskname_suffix));
                    outdirname = "{0}_{1:02d}{2:02d}{3}".format(taskname_prefix, cen1_dst, cen2_dst, taskname_suffix);
                else:
                    outfilename = filename.replace(".root", "_{0:02d}{1:02d}_occupancy{2:d}_{3:d}{4}.root".format(cen1_dst, cen2_dst, occ1, occ2, taskname_suffix));
                    outdirname = "{0}_{1:02d}{2:02d}_occupancy{3:d}_{4:d}{5}".format(taskname_prefix, cen1_dst, cen2_dst, occ1, occ2, taskname_suffix);

                for icen_src in range(0, ncen_src):
                    cen1_src = list_centmin_src[icen_src];
                    cen2_src = list_centmax_src[icen_src];

                    if (cen1_dst <= cen1_src and cen1_src <= cen2_dst) and (cen1_dst <= cen2_src and cen2_src <= cen2_dst):
                        if occ1 == "" and occ2 == "":
                            taskname = "{0}_{1:02d}{2:02d}{3}".format(taskname_prefix, cen1_src, cen2_src, taskname_suffix);
                            list_original_tasknames.append(taskname);
                        else:
                            taskname = "{0}_{1:02d}{2:02d}_occupancy{3:d}_{4:d}{5}".format(taskname_prefix, cen1_src, cen2_src, occ1, occ2, taskname_suffix);
                            list_original_tasknames.append(taskname);

                print(filename, list_original_tasknames, outfilename, outdirname);
                merge(filename, list_original_tasknames, outfilename, outdirname);

#__________________________________________________________
#__________________________________________________________
if __name__ == "__main__":
    filename = "AnalysisResults_HL_243092.root";
    #taskname_prefix = "dimuon_standalone";
    taskname_prefix = "dimuon_global";
    merge_centrality(filename, taskname_prefix);

    filename_all = filename.replace(".root", "_*.root");
    filename_dst = filename.replace(".root", "_merged_centrality_global.root");
    #filename_dst = filename.replace(".root", "_merged_centrality_standalone.root");
    command_hadd_centrality = "ls {0} {1} | xargs hadd -f {2}".format(filename, filename_all, filename_dst);
    print(command_hadd_centrality);
    subprocess.run(command_hadd_centrality, shell=True, stdout=subprocess.PIPE);





    #merge_occupancy(filename_dst, taskname_prefix);
    #filename_all2 = filename_dst.replace(".root", "_*.root");
    #filename_dst2 = filename_dst.replace(".root", "_occupancy.root");
    #command_hadd_occupancy = "ls {0} {1} | xargs hadd -f {2}".format(filename_dst, filename_all2, filename_dst2);

    #print(command_hadd_occupancy);
    #subprocess.run(command_hadd_occupancy, shell=True, stdout=subprocess.PIPE);
