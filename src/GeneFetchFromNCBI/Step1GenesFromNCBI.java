package GeneFetchFromNCBI;

import AoUtils.AoFile;
import AoUtils.AoMath;
import org.apache.commons.lang3.StringUtils;
import pgl.infra.utils.PStringUtils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.util.*;

public class Step1GenesFromNCBI {
    public Step1GenesFromNCBI(){

        /**
         * 4: 获取小麦已克隆的基因 NCBI blast
         */

//        this.geneDeal();

    }

    public void geneDeal(){

        String inGenebankS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/045_geneClone/01_ori/sequence.txt";
        String locusTitleS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/045_geneClone/002/locusTitleDB.txt";
        String locusTitleS2 = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/045_geneClone/002/004_locusTitleDB_keepTa.txt";
        String geneCategory = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/045_geneClone/002/007_geneCategory.txt";
        String locusTitleS3 = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/045_geneClone/002/004_locusTitleDB_keepTa_unique.txt";
//        this.getGeneListFromNCBI();
//        this.getGeneListFromNCBI2();
//        this.removeDuplicate();
//        this.addCNTtooriFasta();
//        this.getSubsetGeneFasta();
//        this.extractLocusTitle(inGenebankS,locusTitleS);
//        this.extractLocusTitleOrganism(inGenebankS,locusTitleS2);
//        this.mkgeneCategory(locusTitleS2,geneCategory);
        AoMath.countCaseInGroup("/Users/Aoyue/Documents/007_geneCategory.txt",1);

    }

    /**
     *
     * @return
     */
    private String getGeneCatogreybyTitle(String title){
        String out = null;

        String[] key1Array = {"powdery mildew","stem rust","brown rust","stripe rust","yellow rust","leaf rust",
                "fusarium head blight","nucleotide-binding site-leucine-rich repeat resistance","fusarium","fungal pathogen",
                "puccinia","hessian fly","nematode","broad-spectrum","disease resistance","pathogen"};

        // Al 按照拆分字符串来判断
        String[] key2Array = {"salt","host","salinity","aluminum","cold","na-cl","drought","temperature",
                "abiotic stress","heat stress","water-stress","ionic stress","osmotic stress","light-stressed",
                "cupric","salinity","osmotic","waterlogging","water","CBF gene"};

        Arrays.sort(key1Array);
        Arrays.sort(key2Array);

        for (int i = 0; i < key1Array.length; i++) {
            if (title.toLowerCase().contains(key1Array[i])){
                out = "Biotic";
                break;
            }else {
                out = "NA";
            }
        }

        if (out.equals("NA")){
            for (int i = 0; i < key2Array.length; i++) {
                if (title.toLowerCase().contains(key2Array[i])){
                    out = "Abiotic";
                    break;
                }else{
                    out = "NA";
                }
            }
        }

        return out;
    }

    private String getGeneCategorybyGene(String element){
        String out = null;
        String genes[] = {"pm","sr","yr","lr","fhb","nbs-lrr","nbs","lrr","nb-arc-lrr","msp"};
        Arrays.sort(genes);
        for (int i = 0; i < genes.length; i++) {
            String gene = genes[i];
            String genetype1 = "(" + genes[i];
            String genetype2 = genes[i] + "-";
            if (element.toLowerCase().startsWith(gene) || element.toLowerCase().startsWith(genetype1)
                    || element.toLowerCase().contains(genetype2)){
                out = "Biotic";
                System.out.println(element);
                break;
            }else {
                out="NA";
            }
        }
        return out;
    }

    public void mkgeneCategory(String infileS, String outfileS){

        HashMap<String,String> hm = new HashMap<>(AoFile.countFileRowNumber_withHeader(infileS)+1);
        String[] geneCategory = {"Domestication", "Yield", "Flowering", "Biotic", "Abiotic", "Other"};

        HashMap<String,String> hmLocusTitle = new HashMap<>(AoFile.countFileRowNumber_withHeader(infileS)+1);

        try {
            BufferedReader br = AoFile.readFile(infileS);
            BufferedWriter bw = AoFile.writeFile(outfileS);
            bw.write("Locus\tTitle\tGeneCategory");bw.newLine();
            String temp = null;
            List<String> l = new ArrayList<>();
            List<String> titleList = new ArrayList<>();

            int cnt = 0;
            while ((temp = br.readLine()) != null) {
                l = PStringUtils.fastSplit(temp);
                cnt++;
                String locus = l.get(0);
                String title = l.get(1);
                hmLocusTitle.put(locus,title);
                titleList = PStringUtils.fastSplit(title," ");
                String category = null;


                if (title.contains(" Al ")){
                    hm.put(locus,geneCategory[4]);
                    System.out.println(title);
                }
                else if (title.toLowerCase().contains("pmm")){
                    hm.put(locus,geneCategory[4]);
                    System.out.println(title);
                }
                else{
                    // based title
                    category = this.getGeneCatogreybyTitle(title);
                    hm.put(locus,category);
                    if (category.equals("NA")){
                        for (int i = 0; i < titleList.size(); i++) {
                            String element = titleList.get(i);
                            //1
//                    if (element.toLowerCase().contains("btr-") ||element.equalsIgnoreCase("q")){
//                        hm.put(locus,geneCategory[0]);
//                        break;
//                    }
                            //2
                            if (element.toLowerCase().contains("vrn") || element.toLowerCase().contains("vernalization") ||
                                    element.toLowerCase().contains("ppd") || element.toLowerCase().contains("photoperiod")){
                                hm.put(locus,geneCategory[2]);
                                break;
                            } else{
                                //Biotic
                                String geneCa = this.getGeneCategorybyGene(element);
                                hm.put(locus,geneCa);
                            }
                        }
                    }
                }
            }

            StringBuilder sb = new StringBuilder();
            for (Map.Entry<String,String> entry: hm.entrySet()){
                String titleS = hmLocusTitle.get(entry.getKey());
                sb.setLength(0);
                sb.append(entry.getKey()).append("\t").append(titleS).append("\t").append(entry.getValue());
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



    public void srcipt_balst(){
        //nohup blastn -db /data4/home/aoyue/vmap2/daxing/software/blastdb/001_wheat/wheatIWGSCv1.0
        // -query sequence_unique.fasta.txt
        // -num_alignments 1
        // -num_threads 90
        // -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen"
        // -out sequence_unique.txt
        // 2>./blast.log &

//        String db = "";
//        String query = "";
//        String num_alignments = "";
//        String num_threads = "";
//        String outfmt = "";
//        String out = "";
//        String log = "";
        String evalue = "";
        String reward = "";

        String db = "";
        String query = "";
        String num_alignments = "";
        String num_threads = "";
        String outfmt = "";
        String out = "";
        String log = "";

        String cmd = "nohup blastn -db " + db + " -query " + query + " -num_alignments " + num_alignments
                + " -num_threads " + num_threads + " -outfmt " + outfmt + " -out " + out
                + " " + log;

        System.out.println(cmd);

    }

    public void extractLocusTitleOrganism(String infileS, String outfileS){

        try {
            BufferedReader br = AoFile.readFile(infileS);
            BufferedWriter bw = AoFile.writeFile(outfileS);
            String temp = null;
            List<String> l = new ArrayList<>();
            int cnt = 0;
            StringBuilder sb = new StringBuilder();
            while ((temp = br.readLine()) != null) {

                if (temp.startsWith("LOCUS")){
                    cnt++;
                    l = Arrays.asList(StringUtils.split(temp," "));
                    String locus = l.get(1);
                    while((temp = br.readLine()) != null){
                        if (temp.startsWith("//")) {
//                            cnt++;
                            break; //意思是我要循环这一个locus 内的东西
                        }

                        if (temp.startsWith("  ORGANISM  ")){ //每一行都要经过 物种的判断,如果不是小麦，就过滤掉。
                            if (!temp.split("  ORGANISM  ")[1].equals("Triticum aestivum"))break;
                        }

                        if (temp.startsWith("  TITLE")){
                            l = PStringUtils.fastSplit(temp,"     ");
                            sb.append(l.get(1)).append(" ");
                            while ((temp = br.readLine()) != null){
                                if (temp.startsWith("  JOURNAL")) break;
                                sb.append(temp.substring(12)).append(" ");
                            }
                            sb.deleteCharAt(sb.length()-1);
                            bw.write(locus + "\t"); bw.write(sb.toString()); bw.newLine();
                            sb.setLength(0);
                            break; //这条语句可以去除第二个或者第3个title， 直接跳出locus循环，直接读下一个locus
                        } // title 后的循环

                    } // locus 后的循环

                }

            }
            br.close();
            bw.flush();
            bw.close();
            System.out.println( cnt + " genes in this db");
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }

    public void extractLocusTitle(String infileS, String outfileS){

        try {
            BufferedReader br = AoFile.readFile(infileS);
            BufferedWriter bw = AoFile.writeFile(outfileS);
            String temp = null;
            List<String> l = new ArrayList<>();
            int cnt = 0;
            StringBuilder sb = new StringBuilder();
            while ((temp = br.readLine()) != null) {

                if (temp.startsWith("LOCUS")){
                    cnt++;

                    l = Arrays.asList(StringUtils.split(temp," "));
                    String locus = l.get(1);
                    while((temp = br.readLine()) != null){
                        if (temp.startsWith("//")) {
//                            cnt++;
                            break; //意思是我要循环这一个locus 内的东西
                        }
                        if (temp.startsWith("  TITLE")){
                            l = PStringUtils.fastSplit(temp,"     ");
                            sb.append(l.get(1)).append(" ");
                            while ((temp = br.readLine()) != null){
                                if (temp.startsWith("  JOURNAL")) break;
                                sb.append(temp.substring(12)).append(" ");
                            }
                            sb.deleteCharAt(sb.length()-1);
                            bw.write(locus + "\t"); bw.write(sb.toString()); bw.newLine();
                            sb.setLength(0);
                            break; //这条语句可以去除第二个或者第3个title， 直接跳出locus循环，直接读下一个locus
                        } // title 后的循环

                    } // locus 后的循环

                }

            }
            br.close();
            bw.flush();
            bw.close();
            System.out.println( cnt + " genes in this db");
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }


    public void getSubsetGeneFasta(){
//        String genedbS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/045_geneClone/002/GeneID_removeDuplicate.txt";
        String genedbS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/045_geneClone/002/006_GeneID_removeDuplicate.txt";
        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/045_geneClone/002/003_sequence_addCnt.fasta.txt";
        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/045_geneClone/003/sequence_unique.fasta.txt";

        List<String> geneNameList = new ArrayList<>();
        List<String> fastaList = new ArrayList<>();
        HashMap<String,String> hmIndexGene = new HashMap<>();
        List<String> indexList = new ArrayList<>();

        List<String> locusList = AoFile.getStringListwithoutHeader("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/045_geneClone/002/004_locusTitleDB_keepTa.txt",0);

        try {
            BufferedReader br = AoFile.readFile(infileS);
            String temp = null;
            List<String> l = new ArrayList<>();
            int cnt = 0;
            StringBuilder sb = new StringBuilder();
            while ((temp = br.readLine()) != null) {
                sb.setLength(0);
                if (temp.startsWith(">")){
                    l = PStringUtils.fastSplit(temp," ");
                    String name = l.get(0).replaceFirst(">","");
                    String locus = l.get(1).split("\\.")[0];
                    while((temp = br.readLine()) != null){
                        if (temp.isEmpty() || temp == "" || temp == null) break; //由于最后一行是空
                        sb.append(temp).append("\n");
                    }

                    String fasta = sb.toString();
                    /**
                     * filter fasta length more than 30 kb;;;;;; also filter genes not in triticum
                     */
                    if (fasta.length() > 30000) continue;
                    if (Collections.binarySearch(locusList,locus) < 0) continue; // 说明该基因不在小麦物种里


                    geneNameList.add(name);
                    fastaList.add(fasta);
                    cnt++;
//                    System.out.println(name);
                }
            }
            br.close();

            System.out.println(cnt + " genes ************************************************************************ in the inatial db");


        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }

        try {
            BufferedReader br = AoFile.readFile(genedbS);
            BufferedWriter bw = AoFile.writeFile(outfileS);
            String temp = null;
            List<String> l = new ArrayList<>();
            while ((temp = br.readLine()) != null) {
                l = PStringUtils.fastSplit(temp);
                String indexID = l.get(1);
                String gene = l.get(0) + ":" + l.get(1)+ ":" + l.get(2);
                indexList.add(indexID);
                hmIndexGene.put(indexID,gene);
            }
            System.out.println("total " + indexList.size() + " unique genes finally");

            Collections.sort(indexList);
            int cnt = 0;
            for (int i = 0; i < geneNameList.size(); i++) {
                int index = Collections.binarySearch(indexList,geneNameList.get(i)); //在库中是第几个 index
                if (index > -1){
                    cnt++;
                    bw.write(">" + hmIndexGene.get(indexList.get(index))); bw.newLine();
                    bw.write(fastaList.get(i)); bw.newLine();
                }
            }

            System.out.println(cnt + " genes remain in the final fasta beacuse remove duplicates and 30 kb longer bases");
            bw.flush();
            bw.close();

        }catch (Exception e){
            e.printStackTrace();
            System.exit(1);
        }

    }


    public void removeDuplicate(){
//        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/045_geneClone/002/GeneID.txt";
//        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/045_geneClone/002/GeneID_removeDuplicate.txt";

        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/045_geneClone/002/005_GeneID_addIndex_addLOCUS.txt";
        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/045_geneClone/002/006_GeneID_removeDuplicate.txt";

        try {
            BufferedReader br = AoFile.readFile(infileS);
            BufferedWriter bw = AoFile.writeFile(outfileS);
//            String header = br.readLine();
//            bw.write(header);bw.newLine();
            String temp = null;
            List<String> l = new ArrayList<>();
            Set<String> geneSet = new HashSet<>();
            int cnt = 0;
            while ((temp = br.readLine()) != null) {
                l = PStringUtils.fastSplit(temp);
                String gene = l.get(2); //******************** need to change
                if (geneSet.add(gene)){
                    bw.write(temp);
                    bw.newLine();
                }
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

    public void addCNTtooriFasta(){
        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/045_geneClone/01_ori/sequence.fasta.txt";
        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/045_geneClone/002/sequence_addCnt.fasta.txt";
        try {
            BufferedReader br = AoFile.readFile(infileS);
            BufferedWriter bw = AoFile.writeFile(outfileS);
//            String header = br.readLine();
//            bw.write(header);bw.newLine();
            String temp = null;
            List<String> l = new ArrayList<>();
            int cnt = 0;
            StringBuilder sb = new StringBuilder();
            while ((temp = br.readLine()) != null) {
                sb.setLength(0);
                if (temp.startsWith(">")){
                    cnt++;
                    sb.append(">Index").append(String.valueOf(cnt)).append(" ");
                    String head = sb.toString();
                    bw.write(head);bw.write(temp.substring(1));
                    bw.newLine();
                }
                else{
                    bw.write(temp);
                    bw.newLine();
                }
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
     * extract gene name and find unique one for blast
     */
    public void getGeneListFromNCBI2(){
        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/045_geneClone/01_ori/GeneFullID.txt";
        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/045_geneClone/002/GeneID_addIndex_addLOCUS.txt";
        try {
            BufferedReader br = AoFile.readFile(infileS);
            BufferedWriter bw = AoFile.writeFile(outfileS);
//            String header = br.readLine();
//            bw.write(header);bw.newLine();
            String temp = null;
            List<String> l = new ArrayList<>();
            int cnt = 0;
            while ((temp = br.readLine()) != null) {
                l = PStringUtils.fastSplit(temp);
                cnt++;
                if (temp.contains("(")){ // cycle1
                    int start = temp.indexOf("(");
                    int end = temp.indexOf(")");
                    String sub1 = temp.substring(start+1,end);

                    String locus = temp.substring(0,20).split("\\.")[0].replace(">","");
//                    System.out.println(sub1);
                    bw.write(locus + "\t" + "Index"); bw.write(String.valueOf(cnt));bw.write("\t");bw.write(sub1);
                    bw.newLine();
                    if (temp.contains("allele")){
                        if (temp.contains("gene,")){
                            int start2 = temp.indexOf("gene,");
                            int end2 = temp.indexOf("allele");
                            String sub2 = temp.substring(start2+5,end2);
                            System.out.println(sub2);
//                            bw.write(sub2);bw.write("\t");
                        }
                        if (temp.contains("mRNA,")){
                            int start2 = temp.indexOf("mRNA,");
                            int end2 = temp.indexOf("allele");
                            String sub3 = temp.substring(start2+5,end2);
                            System.out.println(sub3);
//                            bw.write(sub3);bw.write("\t");
                        }
                    }
                } // cycle1
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
     * extract gene name and find unique one for blast
     */
    public void getGeneListFromNCBI(){
        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/045_geneClone/01_ori/GeneFullID.txt";
        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/045_geneClone/002/GeneID.txt";
        try {
            BufferedReader br = AoFile.readFile(infileS);
            BufferedWriter bw = AoFile.writeFile(outfileS);
//            String header = br.readLine();
//            bw.write(header);bw.newLine();
            String temp = null;
            List<String> l = new ArrayList<>();
            int cnt = 0;
            while ((temp = br.readLine()) != null) {
                l = PStringUtils.fastSplit(temp);
                cnt++;
                if (temp.contains("(")){
                    int start = temp.indexOf("(");
                    int end = temp.indexOf(")");
                    String sub1 = temp.substring(start+1,end);
//                    System.out.println(sub1);
                    bw.write("Index ");bw.write(String.valueOf(cnt));bw.write("\t");bw.write(sub1);
                    bw.newLine();
                    if (temp.contains("allele")){
                        if (temp.contains("gene,")){
                            int start2 = temp.indexOf("gene,");
                            int end2 = temp.indexOf("allele");
                            String sub2 = temp.substring(start2+5,end2);
                            System.out.println(sub2);
//                            bw.write(sub2);bw.write("\t");
                        }
                        if (temp.contains("mRNA,")){
                            int start2 = temp.indexOf("mRNA,");
                            int end2 = temp.indexOf("allele");
                            String sub3 = temp.substring(start2+5,end2);
                            System.out.println(sub3);
//                            bw.write(sub3);bw.write("\t");
                        }
                    }
                }
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
