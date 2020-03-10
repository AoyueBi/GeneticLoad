package Plot;

import AoUtils.AoFile;
import pgl.utils.PStringUtils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

public class GOanalysis {

    public GOanalysis(){
//        this.changeGeneVersion();
        this.mkGO2gene();
//        this.addOnt();
//        this.getCytoscapeInput();

    }

    public void getCytoscapeInput(){
        try {
            String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/104_XPCLR/007_GO/001_input/TERM2GENE_v1__HCgenes_v1.0_repr.TEcleaned.txt";
            String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/104_XPCLR/007_GO/001_input/GENE2TERM_v1__HCgenes_v1.0_repr.TEcleaned_cytoscapeInput.txt";
            BufferedReader br = new AoFile().readFile(infileS);
            BufferedWriter bw = new AoFile().writeFile(outfileS);
            bw.write("(species=sss_aa)(type=consortium)(curator=Go)");
            bw.newLine();
            String temp = null;
            List<String> l = new ArrayList<>();
            int cnt = 0;
            while ((temp = br.readLine()) != null) {
                l = PStringUtils.fastSplit(temp);
                cnt++;
                String go = l.get(0).split(":")[1];
                String trans = l.get(1);
                bw.write(trans + " = " + go);
                bw.newLine();
            }
            br.close();
            bw.flush();
            bw.close();
            System.out.println();
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }

    /**
     * 向结果文件中添加Gene Ontology分类信息
     */
    public void addOnt(){
        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/104_XPCLR/007_GO/002_out/test.txt";
        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/104_XPCLR/007_GO/002_out/002_test_addOnt.txt";
        String ontfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/104_XPCLR/007_GO/001_input/TERM2ont_v1__HCgenes_v1.0_repr.TEcleaned.txt";
        HashMap<String,String> hm = AoFile.getHashMapStringKey_withoutHeader(ontfileS,0,1);

        AoFile.addColumbyString(infileS,1,hm,"ONTOLOGY");

    }

    /**
     *
     */
    public void mkGO2gene(){
        try {
            String infileS = "/Users/Aoyue/Documents/Data/wheat/gene/iwgsc_refseqv1.0_FunctionalAnnotation_v1/iwgsc_refseqv1.0_FunctionalAnnotation_v1__HCgenes_v1.0-repr.TEcleaned.TAB";
            String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/104_XPCLR/007_GO/001_input/TERM2gene_v1__HCgenes_v1.0_repr.TEcleaned.txt";
            String outfileS2 = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/104_XPCLR/007_GO/001_input/TERM2name_v1__HCgenes_v1.0_repr.TEcleaned.txt";
            String outfileS3 = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/104_XPCLR/007_GO/001_input/TERM2ont_v1__HCgenes_v1.0_repr.TEcleaned.txt";


            BufferedReader br = new AoFile().readFile(infileS);
            BufferedWriter bw = new AoFile().writeFile(outfileS);
            BufferedWriter bw2 = new AoFile().writeFile(outfileS2);
            BufferedWriter bw3 = new AoFile().writeFile(outfileS3);

            AoFile.readheader(infileS);
            String temp = br.readLine();
            String GOid = null;
            String name = null;
            String ont = null;
            List<String> l = new ArrayList<>();
            List<String> go = new ArrayList<>();
            List<String> gof = new ArrayList<>();
            int cnt = 0;
            while ((temp = br.readLine()) != null) {
                l = PStringUtils.fastSplit(temp);
                cnt++;
//                String trans = l.get(0).split("\\.")[0];
                String trans = l.get(0);
                String GOl = l.get(7);
                if (GOl.isEmpty() || GOl == null || GOl == "")continue;
                go = PStringUtils.fastSplit(GOl,";");
                //GO:0005739 CC: mitochondrion;GO:0032780 BP: negative regulation of ATPase activity;GO:0042030 MF: ATPase inhibitor activity
//GO:0032780 BP: negative regulation of ATPase activity
//                if (go.size() <0.5)continue;
                for (int i = 0; i < go.size(); i++) {
                    GOid = go.get(i).split(" ")[0];
                    name = go.get(i).substring(15);
                    ont = go.get(i).split(" ")[1].substring(0,2);
                    bw.write(GOid + "\t" + trans);
                    bw2.write(GOid + "\t" + name);
                    bw3.write(GOid+"\t"+ont);

                    bw.newLine();
                    bw2.newLine();
                    bw3.newLine();
                }
            }
            br.close();
            bw.flush();
            bw.close();
            bw2.flush();
            bw2.close();
            bw3.flush();
            bw3.close();
            System.out.println();
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }

    /**
     * 将V2版本的基因转换成V1版本的
     */
    public void changeGeneVersion(){
        try {
            String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/104_XPCLR/007_GO/001_input/GeneID.txt";
            String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/104_XPCLR/007_GO/001_input/002_GeneID_V1.txt";
            BufferedReader br = new AoFile().readFile(infileS);
            BufferedWriter bw = new AoFile().writeFile(outfileS);
            String temp = null;
            int cnt = 0;
            while ((temp = br.readLine()) != null) { //19
                cnt++;
                Character goal = temp.charAt(10);
                System.out.println(goal);
                StringBuilder sb = new StringBuilder();
                for (int i = 0; i < temp.length(); i++) {
                    if (i==10){
                        goal = '1';
                        sb.append(goal);
                    }
                    else if (!(i==10)){
                        goal = temp.charAt(i);
                        sb.append(goal);
                    }
                }
                bw.write(sb.toString());
                bw.newLine();
            }
            br.close();
            bw.flush();
            bw.close();
            System.out.println();
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }

}
