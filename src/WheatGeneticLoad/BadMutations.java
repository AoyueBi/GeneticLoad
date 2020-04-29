package WheatGeneticLoad;

import AoUtils.AoFile;
import AoUtils.SplitScript;
import AoUtils.Triads.Triadsgenes;
import pgl.infra.utils.IOFileFormat;
import pgl.infra.utils.IOUtils;
import pgl.infra.utils.PStringUtils;
import pgl.infra.table.RowTable;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

public class BadMutations {

    public BadMutations(){
//        this.getGeneNameList();
//        this.script();
//        this.limitThreads();
//        this.checkIfdone();
//        this.checkIfdone_Prediction();
//        this.getTriadGenesScript();
        SplitScript.splitScript2("/Users/Aoyue/PycharmProjects/SelfFile/vmap2/003_predictionCheck/003_predictionlist_forBadMutations20200428_remained_001_triadGenes_script.sh",10,1157); //11567 genes for triad

    }

    /**
     * 未跑完的 19769 genes 中，先跑是triads gene的那些
     */
    public void getTriadGenesScript(){
        String infileS = "/Users/Aoyue/PycharmProjects/SelfFile/vmap2/003_predictionCheck/002_predictionlist_forBadMutations20200428_remained_script.sh";
        String outfileS1 ="/Users/Aoyue/PycharmProjects/SelfFile/vmap2/003_predictionCheck/003_predictionlist_forBadMutations20200428_remained_001_triadGenes_script.sh";
        String outfileS2 ="/Users/Aoyue/PycharmProjects/SelfFile/vmap2/003_predictionCheck/003_predictionlist_forBadMutations20200428_remained_002_NoneTriadGenes_script.sh";

        Triadsgenes trg = new Triadsgenes();

        try{
            BufferedReader br = AoFile.readFile(infileS);
            BufferedWriter bw = AoFile.writeFile(outfileS1);
            BufferedWriter bw2 = AoFile.writeFile(outfileS2);
            List<String> l = new ArrayList<>();
            String temp = null;
            int cnt =0;
            while((temp = br.readLine()) != null){
                l = PStringUtils.fastSplit(temp," ");
                String gene = l.get(7); ///data4/home/aoyue/vmap2/daxing/analysis/017_badMutation/002_geneFa/TraesCS4D02G288600.1.fasta
                gene = gene.substring(gene.indexOf("Traes"));
                gene = gene.split("\\.")[0];
                boolean out = trg.ifTriads(gene);
                if (out){ // 说明是 triads基因
                    bw.write(temp);
                    bw.newLine();
                    cnt++;
                }else {
                    bw2.write(temp);
                    bw2.newLine();
                }
            }
            br.close();
            bw.flush();
            bw.close();
            bw2.flush();
            bw2.close();
        }catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }

    /**
     * 检查生成的prediction文件是否存在，不存在那么就重新输出脚本
     */
    public void checkIfdone_Prediction(){

        // ori 83165 genes script
        String infileS1 = "/Users/Aoyue/PycharmProjects/SelfFile/vmap2/001_script/sh_vmap2_predict_badmutations_20200107.sh";
        // finished 63396 genes list
        String infileS2 = "/Users/Aoyue/PycharmProjects/SelfFile/vmap2/003_predictionCheck/001_predictionlist_forBadMutations20200428_finishedGenes.txt";
        String outfileS = "/Users/Aoyue/PycharmProjects/SelfFile/vmap2/003_predictionCheck/002_predictionlist_forBadMutations20200428_remained_script.sh";

        List<String> genesFinishedList = AoFile.getStringListwithoutHeader(infileS2,0);

        try{
            BufferedReader br = AoFile.readFile(infileS1);
            BufferedWriter bw = AoFile.writeFile(outfileS);
            List<String> l = new ArrayList<>();
            String temp = null;
            int cnt =0;
            while((temp = br.readLine()) != null){
                l = PStringUtils.fastSplit(temp," ");
                String gene = l.get(7); ///data4/home/aoyue/vmap2/daxing/analysis/017_badMutation/002_geneFa/TraesCS4D02G288600.1.fasta
                gene = gene.substring(gene.indexOf("Traes"));
                gene = gene.replaceFirst(".fasta","");
                int index = Collections.binarySearch(genesFinishedList,gene);
                if (index>-1)continue;
                bw.write(temp);
                bw.newLine();
                cnt++;
            }
            br.close();
            bw.flush();
            bw.close();
        }catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }



    /**
     * 检查生成的tree 和 maf 文件是否达到 83165 个
     */
    public void checkIfdone(){
        String genelistS = "/Users/Aoyue/PycharmProjects/SelfFile/vmap2/genelist_forBadMutations20200111.txt"; //no header
        String checkedFileS = "/Users/Aoyue/PycharmProjects/SelfFile/vmap2/002_step_aligncheck/001_treelist_forBadMutations20200111.txt";
        String outfileS = "/Users/Aoyue/PycharmProjects/SelfFile/vmap2/002_step_aligncheck/002_align_nofinished.txt";

        List<String> genel = new AoFile().getStringListwithoutHeader(genelistS,0);
        List<String> genetree = new AoFile().getStringListwithoutHeader(checkedFileS,0);
        List<String> generemain = new ArrayList<>();
        Collections.sort(genetree);
        try {
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            for (int i = 0; i < genel.size(); i++) {
                int index = Collections.binarySearch(genetree,genel.get(i));
                if(index < 0){
                    bw.write(genel.get(i));
                    bw.newLine();
                }
            }
            bw.flush();
            bw.close();
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }

    public void limitThreads(){

        //no change sum:83165
//        String infileS = "";

//        String infileS = "/Users/Aoyue/PycharmProjects/SelfFile/vmap2/001_script/sh_vmap2_align_badmutations_20200107.sh";
//        String infileS = "/Users/Aoyue/PycharmProjects/SelfFile/vmap2/001_script/sh_vmap2_align_badmutations_20200111.sh";

        //predict
        String infileS = "/Users/Aoyue/PycharmProjects/SelfFile/vmap2/001_script/sh_vmap2_predict_badmutations_20200107.sh";

        //align
//        new SplitScript().splitScript(infileS,"align_",100,832);
//        new SplitScript().splitScript(infileS,"align_",122,8);
        //predict
        new SplitScript().splitScript(infileS,"predict_",130,640);
    }

    public void script (){
        //程序运行时，输入输出路径设置
        //输入文件是2个目录，一个目录是fasta 一个目录是 sub
        String fastaDirS = "/data4/home/aoyue/vmap2/daxing/analysis/017_badMutation/002_geneFa";
        String subDirS = "/data4/home/aoyue/vmap2/daxing/analysis/017_badMutation/003_subs";
        String configS = "/data4/home/aoyue/vmap2/analysis/015_annoDB/015_BAD_Mutations/wheat_Config.txt";
        String outfileDirS = "/data4/home/aoyue/vmap2/analysis/015_annoDB/015_BAD_Mutations/Output_Dir";
        String predictDirS = "/data4/home/aoyue/vmap2/analysis/015_annoDB/015_BAD_Mutations/Predictions_Dir";
        String logPredictDirS = "/data4/home/aoyue/vmap2/analysis/015_annoDB/015_BAD_Mutations/log_Dir/predictlog";
        String logAlignmentDirS = "/data4/home/aoyue/vmap2/analysis/015_annoDB/015_BAD_Mutations/log_Dir/alignmentlog";
        String logcompileDirS = "/data4/home/aoyue/vmap2/analysis/015_annoDB/015_BAD_Mutations/log_Dir/compilelog";

        //脚本路径输出
        String script1S = "/Users/Aoyue/PycharmProjects/SelfFile/vmap2/001_srcipt/sh_vmap2_align_badmutations_20200107.sh";
        String script2S = "/Users/Aoyue/PycharmProjects/SelfFile/vmap2/001_srcipt/sh_vmap2_predict_badmutations_20200107.sh";
        String script3S = "/Users/Aoyue/PycharmProjects/SelfFile/vmap2/001_srcipt/sh_vmap2_compile_badmutations_20200107.sh";

//        String script1S = "/Users/Aoyue/PycharmProjects/SelfFile/vmap2/001_script/sh_vmap2_align_badmutations_20200111hh.sh";


        //本地其他文件
        String genelistS = "/Users/Aoyue/PycharmProjects/SelfFile/vmap2/genelist_forBadMutations20200107.txt"; //no header

        //检查没有跑出来的基因
//        String genelistS = "/Users/Aoyue/PycharmProjects/SelfFile/vmap2/002_step_aligncheck/002_align_nofinished.txt";


        //align
//        try {
//            BufferedReader br = IOUtils.getTextReader(genelistS);
//            BufferedWriter bw = IOUtils.getTextWriter(script1S);
//            String temp = null;
//            int cnt = 0;
//            while ((temp = br.readLine()) != null) {
//                cnt++;
////                temp = temp.replaceFirst(".fasta","").replaceFirst("./","");
//                String fastaS = new File(fastaDirS,temp+ ".fasta").getAbsolutePath();
//                String subS = new File(subDirS,temp + ".subs").getAbsolutePath();
//                String logS = new File(logAlignmentDirS,temp + "_Alignment.log").getAbsolutePath();
//                StringBuilder sb = new StringBuilder();
//                sb.append("/data1/home/aoyue/biosoftware/BAD_Mutations/BAD_Mutations.py -v DEBUG align -c ").
//                        append(configS).append(" -f ").append(fastaS).append(" -o ").append(outfileDirS).
//                        append(" 2> ").append(logS);
//                bw.write(sb.toString());
//                bw.newLine();
//            }
//            br.close();
//            bw.flush();
//            bw.close();
//            System.out.println(cnt +"\tgenes in align cmd");
//        } catch (Exception e) {
//            e.printStackTrace();
//            System.exit(1);
//        }


        //predict
//        try {
//            BufferedReader br = IOUtils.getTextReader(genelistS);
//            BufferedWriter bw = IOUtils.getTextWriter(script2S);
//            String temp = null;
//            int cnt = 0;
//            while ((temp = br.readLine()) != null) {
//                cnt++;
//                temp = temp.replaceFirst(".fasta","").replaceFirst("./","");
//                String fastaS = new File(fastaDirS,temp+ ".fasta").getAbsolutePath();
//                String subS = new File(subDirS,temp + ".subs").getAbsolutePath();
//                String logS = new File(logPredictDirS,temp + "_Predictions.log").getAbsolutePath();
//                String msaS = new File(outfileDirS,temp + "_MSA.fasta").getAbsolutePath();
//                String treeS = new File(outfileDirS,temp + ".tree").getAbsolutePath();
//                StringBuilder sb = new StringBuilder();
//                sb.append("/data1/home/aoyue/biosoftware/BAD_Mutations/BAD_Mutations.py -v DEBUG predict -c ").
//                        append(configS).append(" -f ").append(fastaS).append(" -o ").append(predictDirS).
//                        append(" -a ").append(msaS).append(" -r ").append(treeS).append(" -s ").append(subS).
//                        append(" 2> ").append(logS);
//                bw.write(sb.toString());
//                bw.newLine();
//            }
//            br.close();
//            bw.flush();
//            bw.close();
//            System.out.println(cnt +"\tgenes in predict cmd");
//        } catch (Exception e) {
//            e.printStackTrace();
//            System.exit(1);
//        }

        //compile in py2 env, only this step in py2 env
//        try {
//            BufferedReader br = IOUtils.getTextReader(genelistS);
//            BufferedWriter bw = IOUtils.getTextWriter(script3S);
//            String logS = new File(logcompileDirS,"Compile.log").getAbsolutePath();
//            StringBuilder sb = new StringBuilder();
//            sb.append("/data1/home/aoyue/biosoftware/BAD_Mutations-1.0/BAD_Mutations.py -v DEBUG compile -p ").
//                    append(predictDirS).append(" 2> ").append(logS);
//            bw.write(sb.toString());
//            bw.newLine();
//
//            br.close();
//            bw.flush();
//            bw.close();
//        } catch (Exception e) {
//            e.printStackTrace();
//            System.exit(1);
//        }
    }

    public void getGeneNameList(){
//        System.out.println("ls *.fasta > /data4/home/aoyue/vmap2/analysis/015_annoDB/015_BAD_Mutations/genelist_forBadMutations20200107.txt");

//        System.out.println("find ./ -name '*.fasta' > /data4/home/aoyue/vmap2/analysis/015_annoDB/015_BAD_Mutations/genelist_forBadMutations20200107.txt &");
//        System.out.println("find ./ -name '*.fasta'|sed 's/.fasta//g'|sed 's/\\.\\///g' | head -n 10");
//        System.out.println("find ./ -name '*.fasta'|sed 's/.fasta//g'|sed 's/\\.\\///g' > genelist_forBadMutations20200111.txt");
//
//        System.out.println("find ./ -name '*.tree'|sed 's/.tree//g'|sed 's/\\.\\///g' > ../treelist_forBadMutations20200111.txt");
//        System.out.println("find ./ -name '*_MSA.fasta'|sed 's/_MSA.fasta//g'|sed 's/\\.\\///g' > ../msalist_forBadMutations20200111.txt");

        System.out.println("find ./ -name '*_Predictions.txt'|sed 's/_Predictions.txt//g'|sed 's/\\.\\///g' > ../predictionlist_forBadMutations20200428.txt");
    }
}
