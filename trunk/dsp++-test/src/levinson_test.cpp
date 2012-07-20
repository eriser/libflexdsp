/*!
 * @file levinson_test.cpp
 *
 * @author Andrzej Ciarkowski <mailto:andrzej.ciarkowski@gmail.com>
 */

#include "levinson_test.h"
#include <dsp++/levinson.h>
#include <dsp++/float.h>

const float x[] = {-0.1494026035070419,0.0810016468167305,2.495556831359863,-3.472383499145508,0.04558731243014336,0.3379714488983154,3.556040763854981,-1.314871549606323,3.78404688835144,-3.084386587142944,5.018832206726074,-9.288578033447266,1.302708268165588,-1.029506921768189,9.610445976257324,7.530250549316406,9.905447959899902,-13.06019592285156,1.900242447853088,6.508235931396484,7.023717403411865,-6.940448760986328,2.180450201034546,3.922278642654419,12.05381011962891,-11.54379749298096,-1.249047040939331,2.190834283828735,-1.988669514656067,-1.263354778289795,11.96612167358398,3.914421558380127,4.648869514465332,-3.330785512924194,-6.659256935119629,5.110325336456299,-9.471199989318848,-3.151104211807251,-2.251407623291016,-16.43711280822754,-11.88524532318115,-11.41867637634277,12.6049690246582,17.22013664245606,-5.643626689910889,4.647700309753418,-16.19404411315918,-5.785148143768311,0.2318349182605743,-3.362344264984131,4.930617809295654,-0.09263123571872711,-12.09934711456299,7.819727897644043,2.357759714126587,6.032336235046387,7.098877429962158,0.4055217206478119,7.174777507781982,-0.8407042622566223,-10.98687267303467,2.54589056968689,14.39997005462647,-8.095022201538086,13.7425422668457,-7.462932109832764,7.014951229095459,7.57051420211792,-19.57714653015137,11.97357845306397,13.78938007354736,-3.46355938911438,7.663932323455811,-16.73236656188965,7.312477588653565,-11.24316692352295,12.8388843536377,1.545591473579407,-11.66268253326416,-10.29715251922607,11.71980857849121,-7.221098899841309,17.65194320678711,-2.142554521560669,9.787291526794434,-9.978960037231445,-6.052635192871094,-4.917722225189209,10.40794277191162,1.505387425422669,13.50559425354004,-6.608733654022217,6.325281620025635,11.42012310028076,-16.02392959594727,-5.656399250030518,0.687798023223877,-25.91687965393066,28.73853874206543,-14.00793647766113,8.799826622009277,-11.09903526306152,-26.63277244567871,12.71606349945068,-7.989556789398193,12.67560768127441,-2.986538648605347,-7.705182075500488,3.365824937820435,-5.505523204803467,-26.6031436920166,7.265079975128174,-7.688900470733643,20.40624618530273,3.721698045730591,6.663389205932617,-9.035960197448731,13.68412303924561,-20.36608695983887,40.26667404174805,1.467981934547424,8.467183113098145,4.19928789138794,-12.73675632476807,-15.91454315185547,11.27905654907227,-2.82566499710083,171.6341400146484,-2.82566499710083,11.27905654907227,-15.91454315185547,-12.73675632476807,4.19928789138794,8.467183113098145,1.467981934547424,40.26667404174805,-20.36608695983887,13.68412303924561,-9.035960197448731,6.663389205932617,3.721698045730591,20.40624618530273,-7.688900470733643,7.265079975128174,-26.6031436920166,-5.505523204803467,3.365824937820435,-7.705182075500488,-2.986538648605347,12.67560768127441,-7.989556789398193,12.71606349945068,-26.63277244567871,-11.09903526306152,8.799826622009277,-14.00793647766113,28.73853874206543,-25.91687965393066,0.687798023223877,-5.656399250030518,-16.02392959594727,11.42012310028076,6.325281620025635,-6.608733654022217,13.50559425354004,1.505387425422669,10.40794277191162,-4.917722225189209,-6.052635192871094,-9.978960037231445,9.787291526794434,-2.142554521560669,17.65194320678711,-7.221098899841309,11.71980857849121,-10.29715251922607,-11.66268253326416,1.545591473579407,12.8388843536377,-11.24316692352295,7.312477588653565,-16.73236656188965,7.663932323455811,-3.46355938911438,13.78938007354736,11.97357845306397,-19.57714653015137,7.57051420211792,7.014951229095459,-7.462932109832764,13.7425422668457,-8.095022201538086,14.39997005462647,2.54589056968689,-10.98687267303467,-0.8407042622566223,7.174777507781982,0.4055217206478119,7.098877429962158,6.032336235046387,2.357759714126587,7.819727897644043,-12.09934711456299,-0.09263123571872711,4.930617809295654,-3.362344264984131,0.2318349182605743,-5.785148143768311,-16.19404411315918,4.647700309753418,-5.643626689910889,17.22013664245606,12.6049690246582,-11.41867637634277,-11.88524532318115,-16.43711280822754,-2.251407623291016,-3.151104211807251,-9.471199989318848,5.110325336456299,-6.659256935119629,-3.330785512924194,4.648869514465332,3.914421558380127,11.96612167358398,-1.263354778289795,-1.988669514656067,2.190834283828735,-1.249047040939331,-11.54379749298096,12.05381011962891,3.922278642654419,2.180450201034546,-6.940448760986328,7.023717403411865,6.508235931396484,1.900242447853088,-13.06019592285156,9.905447959899902,7.530250549316406,9.610445976257324,-1.029506921768189,1.302708268165588,-9.288578033447266,5.018832206726074,-3.084386587142944,3.78404688835144,-1.314871549606323,3.556040763854981,0.3379714488983154,0.04558731243014336,-3.472383499145508,2.495556831359863,0.0810016468167305,-0.1494026035070419};
const float a[] = {1,-2.053357484044227e-013,1.073807709417451e-012,-1.147970607462412e-013,4.983791157542328e-013,-4.345412918382863e-013,1.018740647396044e-012,-9.5723429183181e-013,5.311862061319062e-013,-1.382449710263245e-012,1.145084027598387e-012,-5.715428130770306e-013,-3.563815909046753e-014,-1.840971819433435e-012,8.240075288767912e-013,-1.598277066250375e-012,3.734790254839027e-013,-1.792566095559778e-012,1.413424932650287e-012,-4.596323321948148e-013,6.52367049269742e-013,-1.048938713665848e-012,1.072919530997751e-012,7.52065076881081e-013,3.741451592986778e-013,-2.239430862971403e-012,1.574740338128322e-012,-5.712791351086821e-013,1.24411592139495e-012,3.407829574086918e-013,2.367217533105759e-012,-6.095124405192109e-013,2.307043445171075e-013,-5.226929999935237e-013,8.578693311278585e-013,-1.004196725773454e-012,1.716404796070492e-013,-1.845745778439323e-012,9.2656438077654e-013,-1.818101225126156e-012,2.651212582804874e-013,8.317790900491673e-013,1.159572438069745e-012,-1.938005311785673e-012,4.554134847012392e-013,-2.260858167346669e-012,1.499245172453811e-012,-2.555067268872335e-012,-8.215650382226158e-015,1.923017300953234e-012,7.651657085716579e-013,-1.145084027598387e-012,4.876099524153688e-013,-7.215339437038892e-013,2.604805260375542e-012,-2.119748820916811e-012,1.922462189440921e-012,1.147970607462412e-013,-6.692424392440444e-013,-3.194666753358888e-013,1.036948304999896e-012,1.6404100300349e-012,8.644196469731469e-013,-1.902256130392743e-012,2.599143122949954e-012,-3.155364858287157e-012,6.694644838489694e-013,-7.362999099314038e-013,-1.180631981068103e-012,1.345812350450615e-012,-1.01330055457538e-012,-2.23221441331134e-012,3.332667475319795e-012,-1.589839371263224e-012,8.414380303634061e-013,-2.05924166607474e-012,9.034994974399524e-013,-6.03073146976385e-013,-8.094636072542016e-013,1.164734975134252e-012,7.994716000325752e-013,-2.02171612784241e-012,3.478328736150615e-012,-1.970479335255959e-012,1.341038391444727e-012,-1.483813072411522e-012,7.749356711883593e-014,1.006167371642164e-012,-3.885780586188048e-014,4.065636716177323e-013,1.494804280355311e-012,-9.25592935629993e-013,2.512656749331654e-012,-3.310685059432217e-012,2.023270440076885e-012,-1.136646332611235e-012,-1.354805156950079e-012,9.958700530887654e-013,-1.225630708034942e-012,-4.412026299860372e-013,2.643218977027573e-012,-4.372724404788642e-012,3.019140493165651e-012,-2.450040170742796e-012,6.128431095930864e-014,6.213918268827001e-013,-2.362776641007258e-012,-2.55351295663786e-014,1.005195926495617e-012,-4.718447854656915e-015,3.402833570476105e-012,-2.726374681571997e-012,2.128297538206425e-012,-1.272593141976586e-013,-1.219357947945809e-012,5.936362512670712e-013,5.930811397547586e-013,1.054933917998824e-012,5.266898028821743e-013,-1.885491762720903e-012,3.428590744647408e-012,-1.119881964939395e-012,1.438626995309278e-012,-1.656008663530884e-012,-1.338151811580701e-012,1.63247193540883e-012,-2.908007168400673e-012,-2.963006142525648e-012,4.423017507804161e-012,-3.86779497318912e-012,2.540967436459596e-012,-1.307842723008434e-012,3.719247132494274e-013,-7.857048345272233e-013,-1.558531081968795e-012,1.057265386350537e-012,1.316724507205436e-012,-2.074007632302255e-012,1.992850329202156e-012,-1.337596700068389e-012,2.747024829830025e-012,-1.808830862870536e-013,6.827871601444713e-013,1.350586309456503e-012,-1.94844140821715e-012,-7.10709269213794e-013,1.552757922240744e-012,-1.677769034813537e-012,3.963940287121659e-012,-2.274513910549558e-012,1.624922418841379e-012,5.899725152858082e-013,-1.756372824956998e-012,9.070522111187529e-013,-1.325162202192587e-012,-2.069677762506217e-012,2.257249942516637e-012,-4.652944696204031e-012,2.864264381230441e-012,-1.174838004658341e-012,-1.302069563280384e-012,1.320721310094086e-012,-7.682743330406083e-014,-5.871969577242453e-013,4.873879078104437e-013,-1.409983241273949e-012,2.216671290966588e-012,-1.522976189605174e-012,2.995159675833747e-012,1.149080830487037e-013,3.540501225529624e-013,1.487310274939091e-012,-5.517808432387028e-013,4.463096558993129e-014,1.252664638684564e-012,-2.380429187098798e-012,2.455924352773309e-012,-1.551203610006269e-012,9.987566329527908e-013,-1.101341240428155e-013,5.077049891610841e-013,-6.166178678768119e-013,-2.397526621678026e-012,-8.870681966755001e-013,2.555622380384648e-012,-3.783418023317609e-012,1.297774387953865e-012,-1.856736986383112e-012,1.144861982993461e-012,1.090794121694216e-012,-1.405875416082836e-012,1.767475055203249e-013,1.640687585791056e-012,-2.787048369867762e-012,1.548094985537318e-012,-7.544520563840251e-013,3.145927962577844e-012,-7.650546862691954e-013,1.110667113835007e-012,7.424061365668422e-013,-4.25659507641285e-013,-7.405187574249794e-014,1.459055098962381e-012,-1.059152765492399e-012,1.107780533970981e-012,-3.105848911388875e-012,1.318944953254686e-012,1.463273946455956e-013,-2.120525977034049e-013,-1.367794766338193e-013,1.138644734055561e-012,-4.161115896295087e-013,-4.287126209590042e-013,-3.939071291370055e-012,1.377564728954894e-012,-8.01581023779363e-013,-3.90659726789977e-013,-4.066746939201948e-013,1.921351966416296e-012,-1.110334046927619e-012,7.008837954458613e-013,-4.489741911584133e-013,2.309930025035101e-012,-1.176836406102666e-013,2.97761815204467e-013,-1.23351329150978e-012,1.251665437962402e-012,3.045202978668726e-013,7.522871214860061e-013,5.012656956182582e-013,1.942668248489099e-012,-1.587618925213974e-012,3.18856052672345e-013,-1.575184427338172e-012,1.138200644845711e-012,-1.408873018249324e-012,-5.540012892879531e-014,-1.041611241703322e-012,1.150857187326437e-012,-1.093791723860704e-012,7.94475596421762e-013,-1.072697486392826e-012,1.192934639959731e-012,-2.136957277798501e-012,8.266720641358916e-013,-1.455724429888505e-012,3.719802244006587e-013,-5.642153411145046e-013,1.370459301597293e-012,-1.310951347477385e-012,1.094790924582867e-012,-2.327582571126641e-013,1.331823540340338e-012,-2.366995488500834e-013,-1};
const float e = -1.340213637719957e-10;
const float k[] = {0.5421702265739441,24.07405090332031,-0.5309405326843262,-0.8955219984054565,1.029109001159668,5.724104404449463,-1.204044699668884,-0.2260937690734863,-0.145403191447258,0.5634428262710571,-2.227742671966553,-0.4806896150112152,-0.8275208473205566,0.9697796702384949,-152.6861419677734,-0.9855270385742188,0.8896064162254334,-0.04540691897273064,1.225329875946045,-0.7538387775421143,-1.366348385810852,1.397671937942505,10.92534637451172,-1.277894735336304,1.159735798835754,0.7505678534507752,-1.727526426315308,-0.5276936292648315,-1.039458990097046,-1.476792216300964,-2.72063136100769,2.724761247634888,1.030806422233582,-1.544981360435486,-0.4925901591777802,-0.1091790273785591,-0.4277706742286682,-0.3724576532840729,0.4342067539691925,0.6025733351707459,1.708263635635376,-0.1232666522264481,1.681973338127136,-0.4621567130088806,-2.525862455368042,-0.2070374339818955,-1.379936695098877,1.547559022903442,0.3158368170261383,0.9034579396247864,2.021214485168457,-1.587211608886719,-0.5335932374000549,-0.1856754571199417,0.3024731278419495,0.3179955780506134,-0.8608044981956482,-3.474966526031494,0.955694854259491,-8.958456993103027,-0.9968010783195496,-6.810503959655762,1.097916007041931,1.716730237007141,-0.7601805925369263,-0.05763868242502213,1.37243664264679,0.5907681584358215,1.488462448120117,-1.061983466148377,-1.048125505447388,84.15574645996094,1.036315083503723,0.7922458648681641,-0.7638450860977173,-0.6892961263656616,0.4229184985160828,1.101537227630615,-4.789099216461182,-1.805344104766846,1.491692304611206,-1.133413553237915,-1.06389582157135,-0.111984483897686,-2.178136587142944,0.9792998433113098,12.25214958190918,-0.643949031829834,1.909655570983887,-0.915691077709198,0.4871037900447846,-1.024785876274109,-17.02420806884766,1.018365859985352,1.338384985923767,-0.4099706411361694,0.634595513343811,-1.541774034500122,3.074534893035889,2.300313949584961,-0.6075229048728943,0.1706336140632629,-1.012011289596558,-6.114696979522705,0.9942384958267212,-13.31670761108398,-1.023030877113342,3.180671215057373,0.9080634117126465,1.523760318756104,-1.754252552986145,-1.223477959632874,0.2653111815452576,-0.203640803694725,0.3627357184886932,0.4586856365203857,2.648432731628418,-0.9521628618240356,7.333045959472656,1.063576936721802,0.3190267682075501,0.736029863357544,1.555353045463562,-1.624825358390808,-1.08191192150116,-0.224436566233635,-2.837176084518433,1.357854962348938,1.289344310760498,-0.6794014573097229,3.474529027938843,1.618683815002441,-1.214777708053589,0.3834553360939026,-1.157458186149597,-1.503223776817322,0.2491486966609955,-1.084384083747864,-2.158962965011597,2.676188707351685,0.5888378024101257,1.047130465507507,2.027804851531982,-0.1314550936222076,2.236870527267456,-0.5768750309944153,-1.960354328155518,0.2704661786556244,-0.8821845054626465,0.7490732073783875,3.144990682601929,0.1894971877336502,3.054423570632935,-0.950640857219696,0.5805177092552185,0.5614016056060791,0.01619354635477066,-0.4517040550708771,0.8733270764350891,1.867045164108276,0.9610857367515564,78.98108673095703,-0.9984294772148132,0.1425583809614182,-0.1471388787031174,1.364712834358215,2.024007797241211,-0.2220130115747452,1.29928457736969,-0.1019043996930122,0.6482224464416504,2.78276538848877,-1.114489674568176,-0.2509976327419281,-1.091941356658936,1.296380043029785,-3.207431793212891,-1.594059467315674,2.222878456115723,1.41167151927948,1.094137668609619,-1.488385200500488,-0.1942173540592194,-0.4456585049629211,-0.7270358800888062,0.09163244813680649,0.8535248637199402,5.548745632171631,-0.8870771527290344,1.698465347290039,0.5297410488128662,0.1038571149110794,-0.953886091709137,-5.578725337982178,0.7888849973678589,-0.9960019588470459,45.96431350708008,1.000403642654419,-37.84623718261719,-0.9957180023193359,-1.652036070823669,0.447488009929657,1.023154735565186,-27.01151084899902,-1.018715143203735,-0.1332589089870453,-0.2075161784887314,5.492893218994141,-0.03620830923318863,0.8660181760787964,0.9687833189964294,24.97110176086426,-0.9846089482307434,-0.02924603223800659,0.8981781005859375,0.02100546099245548,-0.6288014054298401,-1.161178231239319,2.286548614501953,3.597107887268066,-1.628077983856201,0.9620592594146729,-0.3412777483463287,3.128640413284302,0.9221435189247131,0.154989019036293,-1.203681468963623,-0.7612408995628357,-1.668895840644836,1.531389117240906,3.900169372558594,-1.414265394210815,0.9051635265350342,0.3520078957080841,-1.389203667640686,-1.952568054199219,1.165956497192383,6.013069152832031,-1.012679219245911,0.2593874037265778,0.138133779168129,0.6413586735725403,-1.299558639526367,-1.95153534412384,0.8349617719650269,-1.234542012214661,-1.393471717834473,-0.2972459197044373,0.4025413393974304,6.420641899108887,-0.4813874959945679,1.042821407318115,-0.1913292109966278,-1};

const double R[] = {-0.577780091782170, 1.44160781329426, -0.421972915671950, 5.50561513144785, -0.421972915671950, 1.44160781329426, -0.577780091782170};
const double A[] = {1, -5.77315972805081e-15, 6.60582699651968e-15,	3.44338737683956e-15, -6.49480469405717e-15, 3.99680288865056e-15, -1,00000000000000};
const double K[] = {2.49508045327004, -1.05160634982999, 12.6509686330085, 1.08197594121543, -3.11606978220238, -1.00000000000000};
const double E = 3.35623256279287e-14;

const double A_short[] = {1,	-13.4326009406974,	-2.68057030622505,	12.6509686330086};
const double K_short[] = {2.49508045327004, -1.05160634982999, 	12.6509686330086};
const double E_short = 50.8401698403905;


void dsp::test::levinson_test::test_levinson()
{
	{
		const size_t L = 7;
		dsp::levinson<double> ld(L);
		double aa[L];
		double kk[L - 1];
		double e = ld(R, aa, kk);
		CPPUNIT_ASSERT(dsp::within_range<double>(0.00001)(e, E));
		CPPUNIT_ASSERT(std::equal(A, A + L, aa, dsp::within_range<double>(0.00001)));
		CPPUNIT_ASSERT(std::equal(K, K + L - 1, kk, dsp::within_range<double>(0.00001)));
	}
	{
		const size_t L = 255;
		dsp::levinson<double> ld(L);
		double aa[L];
		double kk[L - 1];
		double e = ld(x, aa, kk);
		CPPUNIT_ASSERT(dsp::within_range<double>(0.00001)(e, ::e));
		CPPUNIT_ASSERT(std::equal(a, a + L, aa, dsp::within_range<double>(0.00001)));
		CPPUNIT_ASSERT(std::equal(k, k + L - 1, kk, dsp::within_range<double>(0.01)));
	}
}

void dsp::test::levinson_test::test_levinson_complex()
{
	const size_t L = 7;
	dsp::levinson<std::complex<double> > ld(L);
	std::complex<double> aa[L];
	std::complex<double> kk[L - 1];
	std::complex<double> e = ld(R, aa, kk);
	CPPUNIT_ASSERT(dsp::within_range<std::complex<double> >(0.00001)(e, E));
	CPPUNIT_ASSERT(std::equal(A, A + L, aa, dsp::within_range<std::complex<double> >(0.00001)));
	CPPUNIT_ASSERT(std::equal(K, K + L - 1, kk, dsp::within_range<std::complex<double> >(0.00001)));
}

void dsp::test::levinson_test::test_levinson_short()
{
	const size_t L = 7;
	const size_t N = 3;
	dsp::levinson<double> ld(L, N);
	double aa[N + 1];
	double kk[N];
	double e = ld(R, aa, kk);
	CPPUNIT_ASSERT(dsp::within_range<double>(0.00001)(e, E_short));
	CPPUNIT_ASSERT(std::equal(A_short, A_short + N + 1, aa, dsp::within_range<double>(0.00001)));
	CPPUNIT_ASSERT(std::equal(K_short, K_short + N, kk, dsp::within_range<double>(0.00001)));
}


void dsp::test::levinson_test::test_levdown()
{
	const size_t P = 3;
	const double anxt[] = {1, 0.5, 1./3, 0.25, 0.2};
	const double ref[] = {1,	0.468750000000000, 0.277777777777778, 0.156250000000000};
	const double ref_k = 0.2;

	double acur[P + 1] = {0};
	double k;

	dsp::levdown(P, anxt, acur, &k, (double*)0, (double*)0);
	CPPUNIT_ASSERT(std::equal(ref, ref + P + 1, acur, dsp::within_range<double>(0.00001)));
	CPPUNIT_ASSERT(dsp::within_range<double>(0.00001)(k, ref_k));
}

void dsp::test::levinson_test::test_levup()
{
	const size_t P = 4;
	const double acur[] = {1, 0.5, 1./3, 0.25, 0.2};
	const double k = 0.717;
	const double ref[] = {1, 0.6434, 0.512583333333333,	0.489, 0.5585, 0.717};
	const double e_ref = 0.000485911;
	const double e = 0.001;

	double anxt[P + 2] = {0};
	double enxt;
 	dsp::levup(P, acur, anxt, k, &e, &enxt);
 	CPPUNIT_ASSERT(std::equal(ref, ref + P + 2, anxt, dsp::within_range<double>(0.00001)));
 	CPPUNIT_ASSERT(dsp::within_range<double>(0.00001)(e_ref, enxt));
}
