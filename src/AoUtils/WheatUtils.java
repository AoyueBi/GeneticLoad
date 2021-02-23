package AoUtils;

import pgl.infra.utils.PStringUtils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

/**
 * @author AoyueBi
 *
 */
public class WheatUtils {

    private static String positionFileS =
                    "1	0	471304005	chr1A	0	471304005\n" +
                    "2	0	122798051	chr1A	471304005	594102056\n" +
                    "3	0	438720154	chr1B	0	438720154\n" +
                    "4	0	251131716	chr1B	438720154	689851870\n" +
                    "5	0	452179604	chr1D	0	452179604\n" +
                    "6	0	43273582	chr1D	452179604	495453186\n" +
                    "7	0	462376173	chr2A	0	462376173\n" +
                    "8	0	318422384	chr2A	462376173	780798557\n" +
                    "9	0	453218924	chr2B	0	453218924\n" +
                    "10	0	348037791	chr2B	453218924	801256715\n" +
                    "11	0	462216879	chr2D	0	462216879\n" +
                    "12	0	189635730	chr2D	462216879	651852609\n" +
                    "13	0	454103970	chr3A	0	454103970\n" +
                    "14	0	296739669	chr3A	454103970	750843639\n" +
                    "15	0	448155269	chr3B	0	448155269\n" +
                    "16	0	382674495	chr3B	448155269	830829764\n" +
                    "17	0	476235359	chr3D	0	476235359\n" +
                    "18	0	139317064	chr3D	476235359	615552423\n" +
                    "19	0	452555092	chr4A	0	452555092\n" +
                    "20	0	292033065	chr4A	452555092	744588157\n" +
                    "21	0	451014251	chr4B	0	451014251\n" +
                    "22	0	222603248	chr4B	451014251	673617499\n" +
                    "23	0	451004620	chr4D	0	451004620\n" +
                    "24	0	58852447	chr4D	451004620	509857067\n" +
                    "25	0	453230519	chr5A	0	453230519\n" +
                    "26	0	256543224	chr5A	453230519	709773743\n" +
                    "27	0	451372872	chr5B	0	451372872\n" +
                    "28	0	261776885	chr5B	451372872	713149757\n" +
                    "29	0	451901030	chr5D	0	451901030\n" +
                    "30	0	114179647	chr5D	451901030	566080677\n" +
                    "31	0	452440856	chr6A	0	452440856\n" +
                    "32	0	165638404	chr6A	452440856	618079260\n" +
                    "33	0	452077197	chr6B	0	452077197\n" +
                    "34	0	268911281	chr6B	452077197	720988478\n" +
                    "35	0	450509124	chr6D	0	450509124\n" +
                    "36	0	23083594	chr6D	450509124	473592718\n" +
                    "37	0	450046986	chr7A	0	450046986\n" +
                    "38	0	286659250	chr7A	450046986	736706236\n" +
                    "39	0	453822637	chr7B	0	453822637\n" +
                    "40	0	296797748	chr7B	453822637	750620385\n" +
                    "41	0	453812268	chr7D	0	453812268\n" +
                    "42	0	184873787	chr7D	453812268	638686055";

    private static String centromereFileS = "chr1A\t210200000\t215800000\n" +
            "chr1B\t237700000\t243500000\n" +
            "chr1D\t166200000\t173800000\n" +
            "chr2A\t326300000\t327000000\n" +
            "chr2B\t344400000\t351300000\n" +
            "chr2D\t264400000\t272500000\n" +
            "chr3A\t316900000\t319900000\n" +
            "chr3B\t345800000\t347000000\n" +
            "chr3D\t237100000\t243200000\n" +
            "chr4A\t264100000\t267900000\n" +
            "chr4B\t303900000\t304400000\n" +
            "chr4D\t182300000\t188200000\n" +
            "chr5A\t252500000\t255100000\n" +
            "chr5B\t198900000\t202500000\n" +
            "chr5D\t185600000\t188700000\n" +
            "chr6A\t283300000\t288700000\n" +
            "chr6B\t323000000\t327500000\n" +
            "chr6D\t211900000\t217400000\n" +
            "chr7A\t360200000\t363800000\n" +
            "chr7B\t308000000\t310100000\n" +
            "chr7D\t336300000\t341700000";


    private static HashMap<String,Integer> hmChromosomeLength = null; //根据染色体号找到该染色体的长度

    private static boolean build = buildMaps ();

    private static boolean buildMaps () {
        String[] temps = positionFileS.split("\n");
        hmChromosomeLength = new HashMap<>();
        String[] temp = null;
        for (int i = 0; i < temps.length; i++) {
            temp = temps[i].split("\t");
            if (i%2 == 0) continue; // i=0,第1行，i除以2没有余数，过滤
            hmChromosomeLength.put(temp[3].replaceFirst("chr", ""), Integer.parseInt(temp[5]));
        }
        return true;
    }


    /**
     * 根据 chr pos 返回标准化到100的pos，取整数。 例如：chr1A 300
     * Note: chr indicate the Ref Chr (1A 1B 1D ...), pos indicate the Ref pos
     * @param chromosome
     * @param posOnchromosome
     * @return
     */
    public static String getScaledPos(String chromosome, int posOnchromosome){
        int chromosomeLength = hmChromosomeLength.get(chromosome);
        return String.format("%.2f",(double)posOnchromosome*100/chromosomeLength);
    }

    /**
     * 根据中心粒的位置，输出一个文件，包含起始终止片段，中间值，以及scale到100后的值
     * @return
     */
    public static File getCentromereFile (String outfileS){
        File out = new File(outfileS);
        String infileS = "/Users/Aoyue/Documents/Data/wheat/position/ChrLenCentPosi_wheat.txt";

        try {
            BufferedReader br = AoFile.readFile(infileS);
            BufferedWriter bw = AoFile.writeFile(outfileS);
            String header = br.readLine();
            bw.write(header + "\tScale100");bw.newLine();
            String temp = null;
            List<String> l = new ArrayList<>();
            int cnt = 0;
            while ((temp = br.readLine()) != null) {
                l = PStringUtils.fastSplit(temp);
                cnt++;
                String chromosome = l.get(0);
                int posOnchromosome = Integer.parseInt(l.get(2));
                String scale = getScaledPos(chromosome,posOnchromosome);
                bw.write(temp + "\t" + scale);
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

        return out;
    }
}
