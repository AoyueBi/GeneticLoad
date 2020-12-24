package Plot;

import AoUtils.AoFile;
import AoUtils.AoMath;

import java.util.HashMap;
import java.util.List;

/**
 * @author AoyueBi
 * @data 2020-09-14 17:30
 */
public class AoMap {
    public AoMap(){
        //        new AoMath().getnlevelsforEachGroup("/Users/Aoyue/Documents/wheatVMapII_germplasmInfo.txt",27,12);
//        this.checkMap();
//        this.mkLRmp();
        this.getAllcountryNum();


    }

    /**
     * 获取整个VMAP2在每个国家的取样数量，画一个整体的分布图
     */
    public void getAllcountryNum(){
        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/001_taxaList/011_taxaInfoDB/taxa_InfoDB.txt";
        List<String> l = AoFile.getStringList(infileS,7);
        AoMath.countCase_fromList_outFile(l);

    }

    /**
     * 查看源Map数据里一共有多少个国家
     */
    public void checkMap(){
        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/map/000_source/map.world.txt";
        String[] region = AoFile.getStringArraybySet(infileS,4);
        System.out.println(region.length + "\tregions in the source data"); //一共252个国家
    }

    /**
     * 做一个只有Landrace的MAP，并
     */
    public void mkLRmp(){
        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/001_taxaList/011_taxaInfoDB/taxa_InfoDB.txt";
        AoFile.readheader(infileS);
        new AoMath().getnlevelsforEachGroup(infileS,15,7);
        HashMap<String,String> hm = new HashMap<>();
//
//        String infileS2 = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/map/002_sourceforMap/003_country_subspecies_20020915.txt";
//        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/map/002_sourceforMap/004_country_subspecies_addLandracebyContinent_20020915.txt";
//        try {
//            BufferedReader br = AoFile.readFile(infileS);
//            String header = br.readLine();
//            String temp = null;
//            List<String> l = new ArrayList<>();
//            int cnt = 0;
//            while ((temp = br.readLine()) != null) {
//                l = PStringUtils.fastSplit(temp);
//
//                String subspecies = l.get(15);
//                String country = l.get(7);
//                String subcontinent = l.get(20);
//                if (!subspecies.equals("Landrace"))continue;
//                hm.put(country,subcontinent);
//                cnt++;
//            }
//            br.close();
//            System.out.println(cnt + "  landraces");
//
//            BufferedWriter bw = AoFile.writeFile(outfileS);
//            br = AoFile.readFile(infileS2);
//            header = br.readLine();
//            bw.write(header + "\tLRbyContinent");bw.newLine();
//            cnt = 0;
//            while ((temp = br.readLine()) != null) {
//                l = PStringUtils.fastSplit(temp);
//                String country = l.get(0);
//                String subcontinent = hm.get(country);
//                cnt++;
//                bw.write(temp + "\t" + subcontinent);
//                bw.newLine();
//            }
//            br.close();
//            System.out.println(cnt + "  countries");
//            bw.flush();
//            bw.close();
//            System.out.println();
//        } catch (Exception e) {
//            e.printStackTrace();
//            System.exit(1);
//        }

    }

}
