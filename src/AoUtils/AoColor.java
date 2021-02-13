package AoUtils;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;

/**
 * @author AoyueBi
 * 
 */
public class AoColor {

    public AoColor(){
        this.getCol();

    }




    /**
     * e.g.
     */
    public void getCol(){
//        String[] in = {"AB","ABD","D"};
//        AoColor.genomeType(in);

//        String[] in = {"Ae.tauschii","Cultivar","Domesticated_emmer","Free_threshing_tetraploid","Landrace","OtherHexaploid","OtherTetraploid","Wild_emmer"};
//        String[] in = {"Cultivar","Domesticated_emmer","Free_threshing_tetraploid","Landrace","OtherHexaploid","OtherTetraploid","Wild_emmer"};
//        String[] in = {"Ae.tauschii","Cultivar","Landrace","OtherHexaploid"};
//        String[] in = {"Domesticated_emmer","Free_threshing_tetraploid","OtherTetraploid","Wild_emmer"};
//        String[] in = {"Ae.tauschii", "Wild_emmer", "Domesticated_emmer","Free_threshing_tetraploid", "OtherTetraploid","Landrace","Cultivar","OtherHexaploid" };
//        String[] in = {"OtherHexaploid","Cultivar","Landrace", "OtherTetraploid","Free_threshing_tetraploid", "Domesticated_emmer","Wild_emmer", "Ae.tauschii"};
//        String[] in = {"OtherHexaploid","Cultivar","Landrace",
//                "OtherTetraploid","Free_threshing_tetraploid","Domesticated_emmer","Wild_emmer","Ae.tauschii"};
//        String[] in = {"Ae.tauschii","Cultivar","Domesticated_emmer","Free_threshing_tetraploid","Landrace","Wild_emmer"};
//        String[] in = {"Ae.tauschii", "Wild_emmer", "Domesticated_emmer","Free_threshing_tetraploid","Landrace","Cultivar",};
        String[] in = { "Wild_emmer", "Domesticated_emmer","Free_threshing_tetraploid","OtherTetraploid","Landrace","Cultivar","OtherHexaploid"};

        AoColor.subspecies(in);

//        String[] in = {"Africa" , "America" ,"Asia" ,   "Europe" , "Oceania"};
//        String[] in = {"Africa", "America", "Central and South Asia", "East Asia", "Europe", "Oceania", "Western Asia"};
//        String[] in = {"Africa", "America", "Central and South Asia", "East Asia", "Europe", "Western Asia"};
//        String[] in = {"America" ,"Europe" , "Western Asia"};
//        String[] in = {"America","Europe","Africa","Western Asia","Central and South Asia","East Asia"};
//        String[] in = {"Europe","America","Western Asia","Africa","Central and South Asia","East Asia"};

//        AoColor.continent_5(in);
//        AoColor.continent_7_shape(in);
//        AoColor.continent_7(in);


//        String[] in = {"2","3","6","5","8","4","1","7","9","15","14","18","12","10","17","11","20","13","16","19"};
//        AoColor.getSubspecies(in);


    }

    public static String getSubspecies(String[] input){
        String out=null;
        HashMap<String,String> hm = new HashMap<>();
//        String[] value = {"Strangulata","Persian wheat","Wild emmer","Domesticated emmer", "Durum",
//                "Ispahanicum","Georgian wheat","Polish wheat","Khorasan wheat","Rivet wheat","Cultivar",
//                "Landrace","OtherHexaploid","Indian dwarf wheat","Macha","Spelt","Tibetan semi-wild","Vavilovii",
//                "Xinjiang wheat","Yunan wheat","Club wheat"};
        String[] value = {"Strangulata (n=36)","Persian wheat (n=9)","Wild emmer (n=66)","Domesticated emmer (n=52)", "Durum (n=18)",
                "Ispahanicum (n=7)","Georgian wheat (n=3)","Polish wheat (n=10)","Khorasan wheat (n=10)","Rivet wheat (n=12)","Cultivar (n=72)",
                "Landrace (n=211)","OtherHexaploid (n=89)","Indian dwarf wheat (n=5)","Macha (n=5)","Spelt (n=14)","Tibetan semi-wild (n=5)","Vavilovii (n=5)",
                "Xinjiang wheat (n=4)","Yunan wheat (n=5)","Club wheat (n=5)"};

        String[] key = new String[21];
        for (int i = 0; i < key.length; i++) {
            key[i] = String.valueOf(i);
            System.out.println(key[i]);
        }
        for (int i = 0; i < value.length; i++) {
            hm.put(key[i],value[i]);
        }

        StringBuilder sb = new StringBuilder();
        sb.append("c(");
        for (int i = 0; i < input.length; i++) {
            String col = hm.get(input[i]);
            sb.append("\"").append(col).append("\",");
        }
        sb.deleteCharAt(sb.length()-1);
        sb.append(")");
        out = sb.toString();
        System.out.println(out);

        return out;
    }

    /**
     *
     * @param input
     * @return
     */

    public static String delVariantsColor(String[] input){
        String out=null;

        HashMap<String,String> hm = new HashMap<>();
        String[] value = {"#d5311c","#e69f00","#004680"};
        String[] key = {"Deleterious","Nonsynonymous","Synonymous"};
        for (int i = 0; i < value.length; i++) {
            hm.put(key[i],value[i]);
        }

        StringBuilder sb = new StringBuilder();
        sb.append("c(");
        for (int i = 0; i < input.length; i++) {
            String col = hm.get(input[i]);
            sb.append("'").append(col).append("',");
        }
        sb.deleteCharAt(sb.length()-1);
        sb.append(")");
        out = sb.toString();
        System.out.println(out);

        return out;
    }

    /**
     *
     * @param input
     * @return
     */

    public static String continent_7_shape (String[] input){
        String out=null;

        HashMap<String,Integer> hm = new HashMap<>();
        Integer[] value = {7,16,8,15,6,18,17};
        String[] key = {"Oceania","Africa","America","Europe","East Asia","Central and South Asia","Western Asia"};
        for (int i = 0; i < value.length; i++) {
            hm.put(key[i],value[i]);
        }

        StringBuilder sb = new StringBuilder();
        sb.append("c(");
        for (int i = 0; i < input.length; i++) {
            int shape = hm.get(input[i]);
            sb.append("").append(shape).append(",");
        }
        sb.deleteCharAt(sb.length()-1);
        sb.append(")");
        out = sb.toString();
        System.out.println(out);

        return out;
    }

    /**
     *
     * @param input
     * @return
     */

    public static String continent_7 (String[] input){
        String out=null;

        HashMap<String,String> hm = new HashMap<>();
        String[] value = {"#7dbde8", "#7B241C", "#dbb3ff","#FF9900","#fc6e6e","#A9A9A9", "#82C782"};
        String[] key = {"Oceania","Africa","America","Europe","East Asia","Central and South Asia","Western Asia"};
        for (int i = 0; i < value.length; i++) {
            hm.put(key[i],value[i]);
        }

        StringBuilder sb = new StringBuilder();
        sb.append("c(");
        for (int i = 0; i < input.length; i++) {
            String col = hm.get(input[i]);
            sb.append("'").append(col).append("',");
        }
        sb.deleteCharAt(sb.length()-1);
        sb.append(")");
        out = sb.toString();
        System.out.println(out);

        return out;
    }


    //分组2为：大洋洲 非洲 北美洲 南美洲 欧洲 亚洲
    //"Oceania","Africa","North America","South America","Europe","Asia"
    //颜色为："#F1E1FF","#F4D03F","#5DADE2","#7B241C","#FF9900","#82C782"
    //
    /**
     *
     * @param input
     * @return
     */

    public static String continent_5 (String[] input){
        String out=null;

        HashMap<String,String> hm = new HashMap<>();
        String[] value = {"#7dbde8", "#7B241C", "#dbb3ff","#FF9900", "#82C782"};
        String[] key = {"Oceania","Africa","America","Europe","Asia"};
        for (int i = 0; i < value.length; i++) {
            hm.put(key[i],value[i]);
        }

        StringBuilder sb = new StringBuilder();
        sb.append("c(");
        for (int i = 0; i < input.length; i++) {
            String col = hm.get(input[i]);
            sb.append("'").append(col).append("',");
        }
        sb.deleteCharAt(sb.length()-1);
        sb.append(")");
        out = sb.toString();
        System.out.println(out);

        return out;
    }

    /**
     *
     * @param input
     * @return
     */

    public static String subspecies (String[] input){
        String out=null;

        HashMap<String,String> hm = new HashMap<>();
        String[] value = {"#87cef9","#9900ff","#7f5701","#016699","#fc6e6e","#fe63c2","#00f3ff","#ffd702"};
        String[] key = {"Ae.tauschii","Cultivar","Domesticated_emmer","Free_threshing_tetraploid","Landrace","OtherHexaploid","OtherTetraploid","Wild_emmer"};
        for (int i = 0; i < value.length; i++) {
            hm.put(key[i],value[i]);
        }

        StringBuilder sb = new StringBuilder();
        sb.append("c(");
        for (int i = 0; i < input.length; i++) {
            String col = hm.get(input[i]);
            sb.append("'").append(col).append("',");
        }
        sb.deleteCharAt(sb.length()-1);
        sb.append(")");
        out = sb.toString();
        System.out.println(out);

        return out;
    }

    /**
     *
     * @param input
     * @return
     */

    public static String subgenome(String[] input){
        String out=null;

        HashMap<String,String> hm = new HashMap<>();
        String[] value = {"#fd8582","#967bce","#4bcdc6"};
        String[] key = {"A","B","D"};
        for (int i = 0; i < value.length; i++) {
            hm.put(key[i],value[i]);
        }

        StringBuilder sb = new StringBuilder();
        sb.append("c(");
        for (int i = 0; i < input.length; i++) {
            String col = hm.get(input[i]);
            sb.append("'").append(col).append("',");
        }
        sb.deleteCharAt(sb.length()-1);
        sb.append(")");
        out = sb.toString();
        System.out.println(out);

        return out;
    }


    /**
     *
     * @param input
     * @return
     */

    public static String genoType2(String[] input){
        String out=null;

        HashMap<String,String> hm = new HashMap<>();
        String[] value = {"#ffd702","#fc6e6e","#87cef9"};
        String[] key = {"AABB","AABBDD","DD"};
        for (int i = 0; i < value.length; i++) {
            hm.put(key[i],value[i]);
        }

        StringBuilder sb = new StringBuilder();
        sb.append("c(");
        for (int i = 0; i < input.length; i++) {
            String col = hm.get(input[i]);
            sb.append("'").append(col).append("',");
        }
        sb.deleteCharAt(sb.length()-1);
        sb.append(")");
        out = sb.toString();
        System.out.println(out);

        return out;
    }

    public static String genomeType(String[] input){
        String out=null;

        HashMap<String,String> hm = new HashMap<>();
        String[] value = {"#ffd702","#fc6e6e","#87cef9"};
        String[] key = {"AB","ABD","D"};
        for (int i = 0; i < value.length; i++) {
            hm.put(key[i],value[i]);
        }

        StringBuilder sb = new StringBuilder();
        sb.append("c(");
        for (int i = 0; i < input.length; i++) {
            String col = hm.get(input[i]);
            sb.append("'").append(col).append("',");
        }
        sb.deleteCharAt(sb.length()-1);
        sb.append(")");
        out = sb.toString();
        System.out.println(out);

        return out;
    }

    public static String models7(String[] input){
        String out=null;

        HashMap<String,String> hm = new HashMap<>();
        String[] value = {"#fd8582","#ffcbcc","#967bce","#cccdfe","darkgray","#4bcdc6","#cdfffc"};
        String[] key = {"M100","M011","M010","M101","M111","M001","M110"};
        for (int i = 0; i < value.length; i++) {
            hm.put(key[i],value[i]);
        }

        StringBuilder sb = new StringBuilder();
        sb.append("c(");
        for (int i = 0; i < input.length; i++) {
            String col = hm.get(input[i]);
            sb.append("'").append(col).append("',");
        }
        sb.deleteCharAt(sb.length()-1);
        sb.append(")");
        out = sb.toString();
        System.out.println(out);

        return out;
    }

    /**
     * 将因子变量进行自然顺序的排序，然后再调用颜色值
     * @param input
     */
    public static void sort(String[] input){
//        String[] order = {"M001","M010","M011", "M100", "M101", "M110", "M111"};
        List<String> l = new ArrayList<>();
        for (int i = 0; i < input.length; i++) {
            l.add(input[i]);
        }
        Collections.sort(l);
        for (int i = 0; i < l.size(); i++) {
            input[i] = l.get(i);
        }
    }

    public static void main(String[] args) {
        new AoColor();
    }

}
