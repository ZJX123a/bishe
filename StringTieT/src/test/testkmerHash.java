package test;

import static org.junit.Assert.*;

import java.util.HashMap;
import java.util.Map;
import java.util.Scanner;
import java.util.Vector;

import org.junit.Before;
import org.junit.Test;

public class testkmerHash {
	kmerHash testKh=new kmerHash();
	Map<Long, Vector<Integer>> testKmerHash = new HashMap<Long, Vector<Integer>>();
	Vector<String> testRead_vector= new Vector<String>();
	@Before
	public void setUp() throws Exception {
	}

	@Test
	public void testReadsToKmer() {
		System.out.println("请输入kmer长度（kmer长度需在15-31之间，默认25）.........");
		Scanner scan = new Scanner(System.in);
		int read = scan.nextInt();
		while (read > 31 || read < 15) {
			System.out.println("您输入的数据不合法，请重新输入！");
			read = scan.nextInt();
		}
		System.out.println("您输入的kmer长度是：" + read);
		testKh.kmer_length = read;
//		testRead_vector.add("CAGCATCAGGTCTCCAGAGCTGCAGAAGACGACGGCCGACTTGGATCACACTCTTGTGAG");
//		testRead_vector.add("TCCGTGAGATCTTCCCAGGGCAGCTCCCCTCTGTGGAATCCAATCTGTCTTCCATCCTGC");
		testRead_vector.add("GGGCCCCTGCCTGGGGGCTTGTCACCTCCCCCACCTTCTTCCTGAGTCATTCCTGCAGCC");
		testRead_vector.add("GAGCCATGCCTAGAGTGGGATGGGCCATTGTTCATCTTCTGGCCCCTGTTGTCTGCATGT");
		testRead_vector.add("CCTTTTCTTTGACCTCTTCTTTCTGTTCATGTGTATTTGCTGTCTCTTAGCCCAGACTTC");
		testRead_vector.add("AGGAGCCATGCCTAGAGTGGGATGGGCCATTGTTCATCTTCTGGCCCCTGTTGTCTGCAT");
		testRead_vector.add("GGGATGGGCCATTGTTCATCTTCTGGCCCCTGTTGTCTGCATGTAACTTAATACCACAAC");
//		testRead_vector.add("CTGGCCATCCGTGAGATCTTCCCAGGGCAGCTCCCCTCTGTGGAATCCAATCTGTCTTCC");
//		testRead_vector.add("AGCATGACTATTTTTAGAGACCCCGTGTCTGTCACTGAAACCTTTTTTGTGGGAGACTAT");
//		testRead_vector.add("CCGTCAGCCTTTTCTTTGACCTCTTCTTTCTGTTCATGTGTATTTGCTGTCTCTTAGCCC");
		testKh.readsToKmer(testKmerHash, testRead_vector);
		for (long key : testKmerHash.keySet()) {
			System.out.println(key+"   "+testKmerHash.get(key));
		}
	/*			377952261816816   [4]
				657438018875516   [4]
				385909140424643   [4]
				245716520777619   [1, 3, 4]
				417736654855948   [4]
				1107543076462331   [1, 3, 4]
				1042737882718046   [0]
				369149081458043   [1, 3, 4]
				982866083110478   [1, 4]
				991619893154747   [2]
				1052489050353119   [2]
				358405752180574   [0]
				1026252294667015   [4]
				744774687162065   [0]
				307129990240139   [0]
				719790330404730   [1, 3, 4]
				984637887183827   [2]
				667036318528852   [0]
				104992500676068   [0]
				713919606532330   [3]
				980574239956701   [2]
				554118380499961   [2]
				545046712581169   [4]
				407019616554941   [1, 3]
				738546962208735   [0]
				832190620757943   [1, 3, 4]
				342904106905060   [1, 3, 4]
				655237223785822   [1, 3, 4]
				1092329903420654   [2]
				1090020509055269   [0]
				1011184009764181   [0]
				921207427692002   [0]
				371076414755799   [0]
				538038050377409   [4]
				178479901633082   [3]
				1052472585321453   [1, 3, 4]
				747586693973495   [0]
				982954525256958   [2]
				793251810344312   [0]
				749846626914685   [0]
				978934442726320   [4]
				534270979151701   [0]
				377783890981981   [0]
				385235657085301   [0]
				415042721498581   [0]
				750411610149983   [0]
				657395926166807   [0]
				984594740101022   [2]
				931204718884809   [2]
				350696418989551   [1, 3, 4]
				553764611914043   [1, 4]
				839448053400340   [4]
				1042310949184840   [2]
				205121820488927   [1, 3]
				1123692687785455   [2]
				375035996810557   [0]
				516020959044307   [0]
				422074513874910   [2]
				832256480884604   [2]
				252650286130494   [1, 3]
				560679239876219   [2]
				347119155011365   [2]
				588779852091116   [2]
				727309458140191   [4]
				276885769115582   [1, 3, 4]
				1105109400811479   [0]
				307723101879673   [0]
				664704857560047   [1, 3]
				1089568490121143   [2]
				1117507405987131   [2]
				527213608024895   [2]
				727298934963013   [0]
				1089158540813548   [4]
				410480216471732   [0]
				344637548243279   [1, 3]
				560851828207438   [2]
				367634363771475   [1, 3]
				102620054117933   [0]
				1046269879136580   [4]
				1105992399916113   [4]
				627361507933673   [1, 3, 4]
				807379796018448   [4]
				279806244114975   [2]
				882814330665939   [1, 3]
				1117071030613949   [2]
				419970002704274   [0]
				1050306852811348   [2]
				153557602135885   [1, 3]
				744639768263253   [1, 3, 4]
				1119224976459901   [2]
				791544076211489   [2]
				977719370388545   [4]
				257646218049445   [1, 3, 4]
				416345460430161   [0]
				603878612444073   [1, 3]
				820487281955709   [1, 3]
				544597239298935   [2]
				742897536022494   [1, 3, 4]
				1090584401927924   [2]
				553980103974473   [0]
				374244080399605   [0]
				726759259367767   [1, 3, 4]
				930150956858233   [1, 3, 4]
				1090573615157223   [2]
				1077226109853170   [2]
				1032027832669463   [0]
				614230408543543   [1, 3]
				823527690717522   [2]
				1054286943482053   [4]
				562398148657019   [2]
				938183929334607   [0]
				557752326913525   [0]
				1010601144521979   [1, 3]
				1116817052662253   [2]
				413278378716879   [2]
				914376491160711   [2]
				502178559377140   [1, 3]
				654858544364180   [1, 3]
				262576713202837   [2]
				1030584872197781   [1, 3, 4]
				1030149314137591   [1, 3]
				539481934878021   [0]
				1076962669346526   [1, 3, 4]
				103319594679219   [2]
				702388035149693   [0]
				*/
		Map<Long, Vector<Integer>> TempMap = new HashMap<Long, Vector<Integer>>();
		Vector<Integer> read_set=new Vector<Integer>();
		read_set.add(1);
		read_set.add(3);
		TempMap.put(163714636091045l, read_set);
		read_set.clear();
		read_set.add(0);
		TempMap.put(702388035149693l, read_set);
	}

	//@Test
	public void testGet_readset() {
		fail("Not yet implemented");
	}

	//@Test
	public void testGet_readset_count() {
		fail("Not yet implemented");
	}

	//@Test
	public void testGet_forward_candidates() {
		fail("Not yet implemented");
	}

	//@Test
	public void testGet_reverse_candidates() {
		fail("Not yet implemented");
	}

	//@Test
	public void testSort_kmer() {
		fail("Not yet implemented");
	}
public void main(String args[]){
	testReadsToKmer();
}
}
