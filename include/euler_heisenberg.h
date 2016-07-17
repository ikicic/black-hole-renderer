#ifndef EULER_HEISENBERG_H
#define EULER_HEISENBERG_H

#include <complex>
#include "numerics.h"

/* REFERENCE: 0107135v2 - QED Effective Action Revisited */

namespace EH {

void generate_lambdas(void);

template <typename _T>
_T __lagrangian_special1(const _T &x) {
  // `x` is proportional to `a`.
  if (x < 0.014) {
    static const real_t coef[9] = {
      -1.0000000000000000000000000000000000,  // x^2
      +6.0000000000000000000000000000000000,  // x^4
      -120.00000000000000000000000000000000,
      +5040.000000000000000000000000000000,
      -362880.00000000000000000000000000000,
      +3.9916800000000000000000000000000000e7,
      -6.227020800000000000000000000000000e9,
      +1.3076743680000000000000000000000000e12,
      -3.5568742809600000000000000000000000e14,
    };
    _T xx = sqr(x);
    return xx * evaluate_polynomial(coef, xx);
  }

  if (x < 0.093) {
    static const real_t num[] = {
      -0.001745944948744094887582964562202,
      -0.1717030664996928155321066,
      -7.216179554658034638373678,
      -172.0291516786873118172506,
      -2594.745364102808504142306,
      -25900.93726248260133538888,
      -172468.1217934906192277872,
      -745728.9974489024390077965,
      -1.933729293345853958157541e6,
      -2.436752192951123055919940e6,
      -569652.3383504729677956362
    };
    static const real_t den[] = {
      1.0000000000000000000000000000000000,
      51.20159426205703943292093,
      1180.126823412350698095126,
      15780.18264386924412802690,
      133179.8153686322214326390,
      725302.7883565152401540210,
      2.504369768276432076839230e6,
      5.150588987497474750523476e6,
      5.491694913879232371951379e6,
      2.188726571126091576493851e6,
      178250.9528761661423908346
    };
    return evaluate_pade(num, den, x - 0.042);
  }

  if (x < 0.25) {
    static const real_t num[] = {
      -0.0207499640022309061904340548393,
      -0.6076274642623898515341822,
      -7.409776820598485255718525,
      -48.54814512149084770841016,
      -182.8340890279561622052278,
      -387.8699797893707233623935,
      -409.6723386513942855809139,
      -146.6651022898971724804342,
      -10.06249955931292974918871
    };
    static const real_t den[] = {
      1.0000000000000000000000000000000000,
      17.27980178492539739735949,
      122.3853049595083657042776,
      451.1445706993403213081461,
      908.2403067763676091285249,
      954.3882086953391594041049,
      453.4773293454778972811257,
      75.91175020849081163488936,
      2.353786482279473262798830
    };
    return evaluate_pade(num, den, x - 0.152);
  }

  if (x < 0.55) {
    static const real_t num[] ={
      -0.09309001442231006289830803506079,
      -1.117143859581950924162603,
      -5.347128278588276378172376,
      -12.68625381907427895591736,
      -14.48936157700478821656792,
      -4.784736453922426720573854,
      3.770626033009502422592029,
      2.324160072160248710751084,
      0.2235728982648905061940671
    };
    static const real_t den[] = {
      1.0000000000000000000000000000000000,
      7.880301523864462315565614,
      23.55124420539295950095869,
      31.77439769298739179080489,
      15.43257986064949894388953,
      -4.582674467753814134197960,
      -6.048295985534199795516016,
      -1.362142458594493049448355,
      -0.05903169386284142847372354
    };
    return evaluate_pade(num, den, x - 0.37);
  }

  if (x < 1.35) {
    static const real_t num[] = {
      -0.285502201761837388033000616468001,
      -1.505297812575724527480732,
      -3.223961317367794382548716,
      -3.582977288949989843277389,
      -2.188124874269145649995882,
      -0.7174207479529588471310896,
      -0.1154980033251606582876626,
      -0.007481981584592237355378942,
      -0.0001171045026151596142362963
    };
    static const real_t den[] = {
      1.0000000000000000000000000000000000,
      3.896650324438007990233407,
      6.083706234837548088507553,
      4.857508473452302720913683,
      2.103206533208570260381563,
      0.4842198470407589125601097,
      0.05402323438240705027533128,
      0.002344414221369615791798811,
      0.00002189835071652344364502848
    };
    return evaluate_pade(num, den, x - 0.85);
  }

  if (x < 3.3) {
    static const real_t num[] = {
      -0.700799943385466779577268083728272,
      -1.44114582223765848781984550,
      -1.20358252248599149877916991,
      -0.52364578836213241599116417,
      -0.12668964988294833934504943,
      -0.016870019651769522781630045,
      -0.0011457429086673159171087058,
      -0.00003306588452004979530692970,
      -2.521173516691485221825485e-7
    };
    static const real_t den[] = {
      1.0000000000000000000000000000000000,
      1.66066361403067384544919281,
      1.11259143309994688878738633,
      0.385763128981767381792016357,
      0.073844812598481078174747857,
      0.0076993831255670107609320303,
      0.00040175423433413469233094405,
      8.5524147134047180678213578e-6,
      4.2705351396928872683425125e-8
    };
    return evaluate_pade(num, den, x - 2.1);
  }

  static const real_t num1[] = {
   +0.00091682535798392418062812800639316598,
   -0.0043336272214974373288276947690818462,
   -0.051806288119930130366922275627197489,
   +0.20005316069810430817817547977755745,
   +0.51363972750982790944355322876571104,
   -1.5830290414348358959897057980788099,
    0.57721566490153286060651209008240243
  };

  static const real_t den1[] = {
    +5.2967474660690024750872200079677279e-6,
    -7.9617003697984854377242420619326793e-6,
    +0.00055181651541820738583069882900436577,
    -0.00065202304568243158159035422542782881,
    +0.032844333483828741064803423663842698,
    -0.021192624150326990579349838848307044,
    1.0000000000000000000000000000000000
  };
  _T result = evaluate_pade(num1, den1, x);

  static const real_t num2[] = {
    -0.00037234226852870920667530837022362446, // x^0
    +0.027388289676425269645608628659476117,  // x^2
    -0.47059578839239856189008731381612738,
    1.0000000000000000000000000000000000
  };

  static const real_t den2[] = {
   3.2355434897807779163711367101197610e-6,  // x^0
   0.00042372881355932203389830508474576271, // x^2
   0.029404211607601438109912686183872625,
   1.0000000000000000000000000000000000
  };

  using std::log;
  result -= log(x) * evaluate_pade(num2, den2, sqr(x));
  return result;
}


template <typename _T>
_T __lagrangian_special2(const _T &x) {
  // `x` is proportional to `b`.
  if (x < 0.014) {
    static const real_t coef[] = {
      2.0000000000000000000000000000000000,  // x^2
      12.000000000000000000000000000000000, // x^4
      240.00000000000000000000000000000000,
      10080.000000000000000000000000000000,
      725760.0000000000000000000000000000,
      7.983360000000000000000000000000000e7,
      1.2454041600000000000000000000000000e10,
      2.6153487360000000000000000000000000e12,
      7.113748561920000000000000000000000e14,
    };
    _T xx = sqr(x);
    return xx * evaluate_polynomial(coef, xx);
  }

  if (x < 0.026) {
    static const real_t coef[] = {
      0.00080193562582911696825838426161444599647770,
      0.0803887151471442471840096540161051584952,
      2.029394958633313448250334844338481068,
      1.000327877398853891665313968762429,
      13.563813646148203672289582497942,
      33.997808109493341251315097025,
      383.252392158226473300595911,
      2589.124637647861956844568,
      33277.816907477098166726,
      439970.9183724416668594,
      8.2202412184962910051e6,
      1.948144802076476452e8,
      5.77391610813747148e9,
      1.6407583064595755e11,
      1.70011360131167e12
    };
    return evaluate_polynomial(coef, x - 0.02);
  }

  if (x < 0.042) {
    static const real_t coef[] = {
      0.0023284265221913011602231801304050943827434,
      0.1379567832624169524760886447995975216281,
      2.0885523448478108884815473857950473396,
      1.8525128119036160296947646215785464,
      17.46254310594872556125859089305610,
      86.31491760948241339474107333418,
      1062.373055445438428482258767118,
      15282.088921289591531121219299,
      210013.2702316641777313511463,
      -1.21667608490607229079094698e6,
      -1.8552595580234317119255128e8,
      -3.810970464138287906501539e9,
      6.45291701007569924490230e10,
      3.5794792785598643236522e12,
      -2.218978002801680444273569e13,
      -2.745006853628250340061775e15,
      1.97641460481620691672e16,
      2.00670941288245571334e18,
      -3.024187659981441135756728e19,
      -1.244758694632048401994767e21,
      3.92657145453468363e22
    };
    return evaluate_polynomial(coef, x - 0.034);
  }

  if (x < 0.068) {
    static const real_t coef[] = {
      0.00616759734679463277576018257607787134681185,
      0.22889379908777364968868099195577654007327,
      2.26394017611781128754997478339385031637,
      3.996449654374204152501555845153996117,
      37.300606017782366417671799849059873,
      262.1573125297992246138613100723261,
      -992.31031245040000786955106516928,
      -79806.520597087330427859432446786,
      -506283.1776640698679963170797460,
      1.826620659637555425988692883911e7,
      1.40016459171261718193461548241e8,
      -5.3616666819327394026358662509e9,
      -7.55249089804591762191516575e9,
      1.561216053022537476446139460e12,
      -1.371171936800677667614098616e13,
      -2.84109971315931504716691315e14,
      7.7817274868402305942648087e15,
      -3.213346965510810290460921e16,
      -1.773214305198916697515641e18,
      3.97393473395245691044291e19,
      -1.883972067137563115536871e20
    };
    return evaluate_polynomial(coef, x - 0.055);
  }

  if (x < 0.105) {
    static const real_t coef[] = {
      0.0155969358902199247020925883771330,
      0.38579446022124298946228638756005,
      2.871640231463446909838971287752,
      8.39767261028470553252662680832,
      3.962404097378832408508720227,
      -738.7193241252103449457413632,
      -2513.785662133512658225701096,
      72500.2383772043602047354896,
      5983.078223779626426786747,
      -7.6610753499188575820641046e6,
      5.6126436000422644225492522e7,
      3.507546096461670074546095e8,
      -9.849095916340582873382294e9,
      6.88635611570746142370799e10,
      3.19419556409654641507945e11,
      -1.212339895904730374672466e13,
      1.24079743397184751014937e14,
      -2.826388266337209986349660e14,
      -1.082776499999958753939627e16,
      1.9569507716690286994419e17,
      -1.682914160321223934878691e18
    };
    return evaluate_polynomial(coef, x - 0.086);
  }

  if (x < 0.154) {
    static const real_t coef[] = {
      0.03806968159812984141297252953902657049,
      0.668458620029003678180412621253584113,
      3.4625814560028314510307040320263388,
      -2.095243756615504569885812012876079,
      -91.11142101222060852105184680985680,
      84.882078519800486122094480305166,
      3444.94319344221925200017240228405,
      -19928.1695568306842798054864375556,
      -34442.14563787846231246130193220,
      1.1757064960159126771415778819671e6,
      -7.683306123291724400071932987419e6,
      8.50999192980951076908458625110e6,
      3.06438287383902660522822192412e8,
      -3.40174511280846771535050990212e9,
      1.85877654178626113296500566123e10,
      -1.7341177931504498120962793582e10,
      -7.9043970553253042560456085936e11,
      9.659389286723624877923852897e12,
      -6.696864419619382370545708157e13,
      2.524649316599230282963043726e14,
      6.61971874584489359174419263e14
    };
    return evaluate_polynomial(coef, x - 0.129);
  }

  if (x < 0.228) {
    static const real_t coef[] = {
      0.0901818317775667696764043461506625626084,
      1.00133840157878013521425318241312481874,
      1.6288516757502872835288330331716746897,
      -14.019104100751003011262606204660315344,
      -2.50664021327063452534447351707664577,
      233.2791643220376652564651934695520245,
      -913.791027496354241324027986352307357,
      -1.31529348909527466338806797934179,
      19016.5888362421758368857899897387395,
      -119419.999579132615139291351157269460,
      379216.40423297705857826941609134982,
      32474.0106380335517323226494919844,
      -8.991778944282352941832995495812890e6,
      6.744732415990998345930847650777005e7,
      -3.123067397874982939752409090024411e8,
      8.72292600292690046484570435583294e8,
      4.7182944680735341057420251019526e8,
      -2.52261847327277270299603321215892e10,
      2.00347453396074380235195422044297e11,
      -1.08060098477703921859089447219979e12,
      4.3145216687027823212606159031276e12
    };
    return evaluate_polynomial(coef, x - 0.190);
  }

  if (x < 0.332) {
    static const real_t coef[] = {
      0.18105511708885602387026194516648988332066,
      1.0024963595092941700767882715848046430620,
      -1.275352302779927180518024117136054911123,
      -6.84383478564445241128502178889546422364,
      25.0651859176868561609489338264642814502,
      -27.686149393887899388552828884124394370,
      -109.453651675056971995591775130196224597,
      729.02079151647162607857144968867435545,
      -2345.0563422901131008133241934332628148,
      4260.150821824199481110233384912253295,
      2326.44187229014717778082733707184183,
      -57953.1112869826828578042025624394230,
      295377.0432973736372138569890983281523,
      -1.033449855230367355377401832370254510e6,
      2.66709154969973753968121404471922134e6,
      -3.9830003177354656295956975303590670e6,
      -6.698275103938557692252378504117053e6,
      8.618534570683219702551530080554058e7,
      -4.4770437807760998883083115925147491e8,
      1.7502814464832855976500607100961439e9,
      -5.613584820331131384401946733089152e9
    };
    return evaluate_polynomial(coef, x - 0.277);
  }

  if (x < 0.484) {
    static const real_t coef[] = {
      0.279700971780030249024952690437808523682447,
      0.50461206124268613051596880229984589316399,
      -2.14449621525761910547252235315111832364271,
      0.6040439696194783815833482157703503407172,
      5.784470892099624896014478023882959528219,
      -18.40553594066592574849682818942977399261,
      33.86516102916696062676289423326968805376,
      -34.0184602066520470307767957421129738561,
      -33.443505789826945444719300416703861774,
      285.079585458142249482045683554867715132,
      -929.65240533769524552556838939290417132,
      2235.6882744925151066009735311307686027,
      -4260.6925844673029440864198116915372870,
      5889.604910975820719517423993577490831,
      -2175.82601730873857150063257023754701,
      -21976.99847098163251486794100516424892,
      104981.1839772170105884019090480213397,
      -332062.681369633534599322992307653209,
      867334.11397739266418687935800058148,
      -1.96804491873159660697772780188876709e6,
      3.8738617877652644702516465764033605e6
    };
    return evaluate_polynomial(coef, x - 0.405);
  }

  if (x < 0.74) {
    static const real_t coef[] = {
      0.3057034633138983277406587753800570277746689,
      -0.180070581689694782250463021315053084531203,
      -1.298248631929604654039896602669135955663007,
      1.50969949681273303842843248308109821905184,
      -0.7818795636068861342912587250765989286768,
      -0.9421754526632936478601636476964401547154,
      3.678197897715607648420551740524844976321,
      -7.215732709150528621938732872500947360875,
      10.90547629814134721360916698778758419493,
      -13.38054416522934727428425023321130492253,
      12.2111530895911806392307008610632916420,
      -3.5286674210567562756744365637988960477,
      -18.283048408922731590747123247505995170,
      60.624736562830203909564102900363215235,
      -132.06994028715386275302622225197793367,
      240.45948697788230974200631850875868318,
      -388.6660976618272226370087788979739805,
      566.4540142573412322343009414878433507,
      -736.172392460184836339767288638446920,
      809.247551847802717431336012982150521,
      -609.76115527167312994169375748469316
    };
    return evaluate_polynomial(coef, x - 0.6);
  }

  if (x < 1.04) {
    static const real_t coef[] = {
      0.1749966395712209497073028424954765168740401,
      -0.6360063792999692469268775710939931278452833,
      -0.408396153726274305550473761904201828269122,
      0.606066194999356984750019799993393529894494,
      -0.540051402549988119140676212644487173951085,
      0.36214123656980838715936044206150975026865,
      -0.14104836768727266293676987123208496316522,
      -0.0837358831654465747097485128003245661565,
      0.2869593399483062675159881528154137893063,
      -0.4519374265924061204084962527364694636936,
      0.567932321496743340168617866051082557258,
      -0.628750265663418812439751699408807777173,
      0.631893933196147048538832380043592151967,
      -0.57799049665870634703865394248766962360,
      0.47036331664237175974664289610240614200,
      -0.31467896001282955264057493754672465913,
      0.1186316709435725474582951350096838270,
      0.1083568376543396482198556956479969259,
      -0.3554356428367529885355417966771079480,
      0.610629067997513802495160830425295531,
      -0.861152149980077109051294721454610233
    };
    return evaluate_polynomial(coef, x - 0.89);
  }

  if (x < 1.3) {
    static const real_t coef[] = {
      -0.0245671705648116590233023363935142623322848,
      -0.7599910581442769760886783970648495174472142,
      -0.0875036799438830413957170317836431046343272,
      0.220920791613837823094405564436358411175498,
      -0.199782742301812636010870641248003284292354,
      0.143457706715668470688233800233173437333834,
      -0.088145018342830283645232081036581686655357,
      0.04447355095212510759904281975411548083872,
      -0.01387795463457513208102877975145943096287,
      -0.00544913025436893659157978174832909790928,
      0.01616780230805484617041332916828520123074,
      -0.0208325073378072196089237912312397469146,
      0.0215539908302761738629009785148098902407,
      -0.0199337213238619968821500626673829628121,
      0.0171113523490672541823313324060097505294,
      -0.0138505293630297269177914732645820029930,
      0.010628792868653428292904600337887650546
    };
    return evaluate_polynomial(coef, x - 1.17);
  }

  if (x < 1.65) {
    static const real_t coef[] = {
      -0.2558000774565338353703628320129531773364768,
      -0.7696946759251928988777721101852639517673887,
      0.0334410914512191763956140374079938345953039,
      0.0733480688948079954063677382784906204424323,
      -0.0696349459502076494219763234672749077761440,
      0.047732401303878262549971381376196880952042,
      -0.028549511750428974945516447414587161426269,
      0.015484596030996816156605459607926485954464,
      -0.007574382131677994706541893542641132642641,
      0.003175161016966284282529374088623633001818,
      -0.000927495990347538178167639831897700348827,
      -0.000098134974931432711489285628140699888261,
      0.000478122477050016788792866486041046451602,
      -0.00054586569719627965798509400577460857509,
      0.00048232810004638902271408824075313355875,
      -0.00037886120948575205878984883969307996178,
      0.00027707011944779963226922526090994677785
    };
    return evaluate_polynomial(coef, x - 1.47);
  }

  if (x < 2.14) {
    static const real_t coef[] = {
      -0.5622235201546276443185463511713115856876453,
      -0.7193046539016755262383386674177810766015597,
      0.0771720314853861845139282265520554158807666,
      0.0114895488337416485354649144301221885663577,
      -0.0176484108052555730586138584700310272544604,
      0.0116977630209470396539638877372210989726749,
      -0.0063914299299164383046963116680563989682885,
      0.0031771388283655872789194270057258974065245,
      -0.0014799762460077290652288030428101583556781,
      0.0006509333147395735694537276606581122602491,
      -0.0002688473545001873584325177702947854507943,
      0.0001020615468163157368722825894336301011685,
      -0.0000336086148300583226823553164009526414475,
      7.7918055784912607329249230818375805485e-6,
      6.284250705589424573952935385890045699e-7,
      -2.5294456078365328264187761173800139552e-6,
      2.3250316259362075177827981589092810952e-6
    };
    return evaluate_polynomial(coef, x - 1.88);
  }

  if (x < 2.8) {
    static const real_t coef[] = {
      -0.94635556466756556688279715620032992960436334,
      -0.6287535988591782072423043373828074037840379,
      0.0769039898983961272429652796825027907446270,
      -0.00627120345133848926463502399286592314969805,
      -0.00229543358266873532161894662762055258202409,
      0.00195121431448123493050904010155177377067053,
      -0.00099752414475625623437841064660342997578440,
      0.00043726642168953481433006141292936013801928,
      -0.00017741223625117712855067552490485295463087,
      0.00006848704054302213227307479182262765526946,
      -0.00002544904472739171014821902984147152424364,
      9.14269671718272788513561579825040537578e-6,
      -3.17494346708582591628918578652400458561e-6,
      1.06102640542146305317096659841968526673e-6,
      -3.37714216100798387327311523597482144586e-7,
      1.00189905816758105720339847007758265492e-7,
      -2.6351190253130675417704900181477429034e-8
    };
    return evaluate_polynomial(coef, x - 2.45);
  }



  static const real_t num1[] = {
    0.00166677657063607666645231665107084368829504106,  // x^0
    0.0954070334489753356963016449053144833504261177,   // x^2
    0.959703231860194451038814287382769779323687739,
    -1.15443132980306572121302418016480486208431867
  };

  static const real_t den1[] = {
    4.52600826574452666644089766454512233609783296e-6,    // x^0
    -0.000512791930298651930399221901307088158181738924,  // x^2
    +0.0319801583763542705226710794131530604648335418,
    -1.00000000000000000000000000000000000000000000
  };

  _T xx = sqr(x);
  _T result = evaluate_pade(num1, den1, xx);

  static const real_t num2[] = {
    0.000744684537057418413350616740447248921825193012, // x^0
    0.0547765793528505392912172573189522342064714946,   // x^2
    0.941191576784797123780174627632254750898818695,
    2.00000000000000000000000000000000000000000000
  };

  static const real_t den2[] = {
    -3.23554348978077791637113671011976096721859434e-6,  // x^0
    0.000423728813559322033898305084745762711864406780,  // x^2
    -0.0294042116076014381099126861838726245505906523,
    +1.00000000000000000000000000000000000000000000
  };

  using std::log;
  result -= log(x) * evaluate_pade(num2, den2, xx);
  return result;
}


template <typename _T, typename _Factor>
_T __lagrangian_special_coth(const _T &x, const _T &y, const _Factor &factor) {
  /* Calculates x * y * coth(factor * x / y) correctly. */
  _T inner = factor * x / y;
  if (factor * x < 0.00001 * y)
    return (sqr(y) / factor) * (1 + sqr(inner) / 3);
  if (factor * x > 1000 * y)
    return x * y;
  // TODO: test coth(x) for autodiff<autodiff<double>> x.
  return x * y * coth(inner);
}

template <typename _T>
_T __lagrangian_exp_inv(const _T &x) {
  /* Calculates exp(-1 / x) for x > 0 correctly. */
  using std::exp;
  if (x < 1e-8)
    return _T(0);  // Not correct in general.
  return exp(-inverse(x));
}

// DO NOT USE std::complex.
// http://stackoverflow.com/q/11108743/2203044
template <typename _T>
std::pair<_T, _T> __lagrangian_ab_helper(const _T &a, const _T &b) {
  using std::abs;
  constexpr double EPSILON = 1e-16;

  _T real = _T();
  _T realc = _T();
  _T imag = _T();
  _T term;
  for (int n = 1; n <= 1000000; ++n) {
    double inv_n_pi = 1 / (M_PI * n);
    _T ax = inv_n_pi * a;
    _T bx = inv_n_pi * b;

    _T cotha = __lagrangian_special_coth(a, b, n * M_PI);
    _T cothb = __lagrangian_special_coth(b, a, n * M_PI);
    // imag += inv_n_pi * cotha * __lagrangian_exp_inv(bx);
    // _T an = cothb * __lagrangian_special1(ax) / n;
    // _T dn = cotha * __lagrangian_special2(bx) * 0.5 / n;
    term = inv_n_pi * (cothb * __lagrangian_special1(ax)
                     - cotha * __lagrangian_special2(bx) * 0.5);
    {
      /* Kahan summation algorithm */
      /* Changes the relative value for 1e-12 when B=2e11T, E=0.  */
#ifdef __FAST_MATH__
#error This method is using Kahan summation algorithm, please make sure \
      that -ffast-math does not interfere with it.
#endif
      _T realy = term - realc;
      _T realt = real + realy;
      realc = (realt - real) - realy;
      real = realt;
    }
    if (n >= 50 && abs(term) < abs(EPSILON * real))
      break;
  }
  // if (abs(term) > 1e-8 * abs(real)) {
  //   std::cerr << "Too low precision: " << term << " "
  //             << real << " " << term / real << '\n';;
  // }

  // using std::isnan;
  // if (isnan(real)) {
  //   std::cerr << "  a=" << a << '\n';
  //   std::cerr << "  b=" << b << '\n';
  //   assert(0 == 1);
  // }
  real *= -1 / M_PI;
  imag *= .5;

  return {real, imag};
}

template <typename _T>
inline std::pair<_T, _T> lagrangian__dimless(_T F, _T G) {
  using std::sqrt;
  F += 1e-20;
  G += 1e-20;
  _T a = sqrt(sqrt(sqr(F) + sqr(G)) + F);
  _T b = G / a;

  std::pair<_T, _T> result = __lagrangian_ab_helper(a, b);
  result.first *= PHY_alpha;
  result.second *= PHY_alpha;
  return result;
}

template <typename _T>
inline std::pair<_T, _T> lagrangian(_T F, _T G) {
  std::pair<_T, _T> result = lagrangian__dimless(
      F / sqr(PHY_Bc),
      G / sqr(PHY_Bc));
  result.first *= sqr(PHY_Bc);
  result.second *= sqr(PHY_Bc);
  return result;
}

template <typename _T>
inline _T lagrangian_real__dimless(const _T &F, const _T &G) {
  return lagrangian__dimless(F, G).first;
}

template <typename _T>
inline _T lagrangian_real(const _T &F, const _T &G) {
  return lagrangian(F, G).first;
}


}

#endif
