#include "tests.h"
#if TESTS_ENABLED
#include "../include/euler_heisenberg.h"
#include "../include/qed_lagrangian.h"

bool test_euler_heisenberg_special1(void) {
  static const double result1[] = {
    0,
    -0.00009994011949958949316934683,
    -0.0003990475545378196175503604,
    -0.0008952243688416267666912836,
    -0.001585101752430308532573549,
    -0.002464206385778246476333385,
    -0.003527153837845227612042549,
    -0.004767846957722380654879948,
    -0.006179663644104410127513577,
    -0.007755624879408931705403825,
    -0.009488539016354807407117484,
    -0.01137112173187010628559611,
    -0.01339609308510601164438021,
    -0.01555625413411956240573804,
    -0.01784454593687884701225455,
    -0.02025409375225369921382930,
    -0.02277823904424392758413961,
    -0.02541056158989846511131823,
    -0.02814489366417184846102455,
    -0.03097532795947663831371313,
    -0.03389622061162176476626653,
    -0.03690219045401092971276601,
    -0.03998811540917151635259778,
    -0.04314912674851048253922930,
    -0.04638060180385888702850387,
    -0.04967815559365675052917451,
    -0.05303763172843484889046913,
    -0.05645509288082107012152250,
    -0.05992681104135016979479534,
    -0.06344925773007957410475267,
    -0.06701909429305281142786066,
    -0.07063316238005060196638523,
    -0.07428847467421951751400659,
    -0.07798220592375553298638196,
    -0.08171168430977622943774436,
    -0.08547438317197499264506695,
    -0.08926791310391344560333334,
    -0.09309001442231006289830804,
    -0.09693855000896882639080103,
    -0.1008114985196962351921222,
    -0.1047069479513841808284153,
    -0.1086230895561543103527758,
    -0.1125582120898774045988894,
    -0.1165106963813474832634763,
    -0.1204790102077835658496914,
    -0.1244617034620553665438406,
    -0.1284574035970053590081684,
    -0.1324648113324070307970068,
    -0.1364826966104088235063280,
    -0.1405098947857264535736663,
    -0.1445453030373324204587029,
    -0.1485878769889264735140644,
    -0.1526366275260358173947717,
    -0.1566906177981742859350805,
    -0.1607489603950743683600326,
    -0.1648108146865862786784349,
    -0.1688753843164078092467292,
    -0.1729419148403627942766877,
    -0.1770096915004812674747060,
    -0.1810780371266485318733703,
    -0.1851463101580819050442784,
    -0.1892139027773620452179485,
    -0.1932802391501901920487489,
    -0.1973447737644634431577728,
    -0.2014069898626576977556727,
    -0.2054663979618827102546566,
    -0.2095225344563265438037185,
    -0.2135749602971384381282678,
    -0.2176232597451106207133177,
    -0.2216670391918118496640599,
    -0.2257059260450994539347905,
    -0.2297395676751933070535433,
    -0.2337676304177354973737344,
    -0.2377897986304843831454230,
    -0.2418057738005021575077114,
    -0.2458152736988918762692142,
    -0.2498180315803239622852469,
    -0.2538137954247642972969698,
    -0.2578023272189769078598251,
    -0.2617834022755246684250576,
    -0.2657568085871320634271131,
    -0.2697223462144055135451628,
    -0.2736798267050296829397824,
    -0.2776290725426731109845364,
    -0.2815699166239439850876985,
    -0.2855022017618373880330006,
    -0.2894257802142093790778716,
    -0.2933405132359012375737329,
    -0.2972462706532195170788383,
    -0.3011429304595546056775822,
    -0.3050303784309926178814463,
    -0.3089085077608429845419903,
    -0.3127772187120673667380785,
    -0.3166364182866547837394585,
    -0.3204860199110433804492580,
    -0.3243259431367413144616313,
    -0.3281561133553480482580667,
    -0.3319764615272231034308391,
    -0.3357869239230922716790221,
    -0.3395874418779215683948066,
    -0.3433779615564270328325330
  };
  if (!evaluate_and_compare(
        [](int i) { return i * 0.01; },
        [](double x) { return EH::__lagrangian_special1(x); },
        result1,
        sizeof(result1) / sizeof(result1[0]),
        2.3e-15)) {
    return false;
  }

  static const double result2[] = {
    -0.3433779615564270328325330,
    -0.3807244450248938272964511,
    -0.4170403262489011549581825,
    -0.4523244714144502446946668,
    -0.4865932786675671280944240,
    -0.5198737649468260301975360,
    -0.5521991076522465591564306,
    -0.5836057503305748747837805,
    -0.6141315202729378756297708,
    -0.6438144091999337261109141,
    -0.6726917928685491115564629,
    -0.7007999433854667795772681,
    -0.7281737376159073644555628,
    -0.7548464971441733785135255,
    -0.7808499162657503459383496,
    -0.8062140484482655652756094,
    -0.8309673310682963139627468,
    -0.8551366345824767160234563,
    -0.8787473266332157402433279,
    -0.9018233445786157609560619,
    -0.9243872720061715639758840,
    -0.9464604162295643679551850,
    -0.9680628847721595472103644,
    -0.9892136595422375651181514,
    -1.009930667894530170839928,
    -1.030230850113251151236195,
    -1.050130223087365438855742,
    -1.069643940110087477511179,
    -1.088786346843224769337779,
    -1.107571033558349293358246,
    -1.126010883811816928295537,
    -1.144118119737162174035965,
    -1.161904344151950642148388,
    -1.179380579680770087648438,
    -1.196557305094543952496529,
    -1.213444489060807093276879,
    -1.230051621491452641160120,
    -1.246387742664779413151171,
    -1.262461470288172466273751,
    -1.278281024656945375050675,
    -1.293854252054111706497489,
    -1.309188646525369594100405,
    -1.324291370153528867184014,
    -1.339169271947077062075532,
    -1.353828905448618833282207,
    -1.368276545160553179894686,
    -1.382518201877574408096662,
    -1.396559637008381931278979,
    -1.410406375962337902139644,
    -1.424063720670691746196970,
    -1.437536761306365298624584,
    -1.450830387261128464196349,
    -1.463949297434260019870693,
    -1.476898009882448953556498,
    -1.489680870876717414006188,
    -1.502302063408507349578502,
    -1.514765615183741425929299,
    -1.527075406140618885628335,
    -1.539235175524114576146621,
    -1.551248528547592201029747,
    -1.563118942669600493851816,
    -1.574849773511774715835600,
    -1.586444260441798463944530,
    -1.597905531843576556476017,
    -1.609236610095114421070327,
    -1.620440416273079896395919,
    -1.631519774601627789269397,
    -1.642477416661785086620854,
    -1.653315985376515565028774,
    -1.664038038785497717852571,
    -1.674646053622651290017972,
    -1.685142428708527869154037,
    -1.695529488168833194831457,
    -1.705809484489566989230626,
    -1.715984601418544608923095,
    -1.726056956722398595762017,
    -1.736028604807542646705882,
    -1.745901539213011421411044,
    -1.755677694982563130048831,
    -1.765358950922944498466612,
    -1.774947131754766307485482,
    -1.784444010162019340812702,
    -1.793851308745872597996221,
    -1.803170701888035611228228,
    -1.812403817528632432587058,
    -1.821552238863224306171603,
    -1.830617505963329354636213,
    -1.839601117324519095774686,
    -1.848504531345921708564913,
    -1.857329167744729265832917,
    -1.866076408909089337011163
  };
  if (!evaluate_and_compare(
        [](int i) { return (i + 10) * 0.1; },
        [](double x) { return EH::__lagrangian_special1(x); },
        result2,
        sizeof(result2) / sizeof(result2[0]),
        6e-16)) {
    return false;
  }

  return true;
}

bool test_euler_heisenberg_special2(void) {
  static const double result1[] = {
    0,
    2.000012000240010080725840e-6,
    8.000192015362581223505443e-6,
    0.00001800097217502617777788877,
    0.00003200307298370136523725770,
    0.00005000750375394460706704025,
    0.00007201556321441458797566568,
    0.00009802884029407531685121243,
    0.0001280492150844591577921169,
    0.0001620788599823043535581687,
    0.0002001202410153387058944283,
    0.0002421761193544572715165983,
    0.0002882495530160409397795404,
    0.0003383438987586886782597202,
    0.0003924628141791936333923359,
    0.0004506102600131872942254840,
    0.0005127905026465124062843669,
    0.0005790081168440709538969602,
    0.0006492679887036360627611176,
    0.0007235753188429251698208414,
    0.0008019356258291169682583847,
    0.0008843547498609692572452389,
    0.0009708388567147744194050913,
    0.001061394441966592936923219,
    0.001156028335504557208630004,
    0.001254747706346568935613374,
    0.001357560067780464210479196,
    0.001464473282845745166860228,
    0.001575495570178346656329670,
    0.001690635510242711594966147,
    0.001809902051978798933852682,
    0.001933304519895665428251482,
    0.002060852621648065574135787,
    0.002192556456138195605373606,
    0.002328426522191301160223185,
    0.002468473727861313714513680,
    0.002612709400430789522239760,
    0.002761145297177854518984493,
    0.002913793616991103587700451,
    0.003070667012920802245346038,
    0.003231778605760507984299916,
    0.003397141998756506481446070,
    0.003566771293542372127852042,
    0.003740681107391699886555361,
    0.003918886591872934415671616,
    0.004101403452975755603287662,
    0.004288247972758424802565156,
    0.004479437032539887046999912,
    0.004674988137629582449660982,
    0.004874919443552443386830277,
    0.005079249783687290555818188,
    0.005287998698194843644720025,
    0.005501186464068035116567944,
    0.005718834126093552579626029,
    0.005940963528470858240211599,
    0.006167597346794632775760184,
    0.006398759120069869336537364,
    0.006634473282396773380941230,
    0.006874765193936110562967581,
    0.007119661170745400719400107,
    0.007369188513062887635122977,
    0.007623375531609817676443486,
    0.007882251571482325770399301,
    0.008145847033212051227755455,
    0.008414193390589209869112273,
    0.008687323204862801332011517,
    0.008965270134959372390511261,
    0.009248068943393629005843617,
    0.009535755497580458285041146,
    0.009828366766297805482506041,
    0.01012594081109254622713871,
    0.01042851677246619462887928,
    0.01073613485072320627984163,
    0.01104883628141101847530566,
    0.01136666330532711333768333,
    0.01168965913311364654860038,
    0.01201786790450397485799434,
    0.01235133464232723040384787,
    0.01269010520141649394009423,
    0.01303422621260275155499067,
    0.01338374502201039282830040,
    0.01373870962590030567280581,
    0.01409916860133349415302312,
    0.01446517103295150224458501,
    0.01483676643618973736119517,
    0.01521400467725607202413806,
    0.01559693589021992470209259,
    0.01598561039156648490937858,
    0.01638007859257698934784059,
    0.01678039090989913861186439,
    0.01718659767467205285542954,
    0.01759874904056780155286309,
    0.01801689489110672067615112,
    0.01844108474659667149260161,
    0.01887136767103732379189772,
    0.01930779217931968814603286,
    0.01975040614503869975423062,
    0.02019925670922288848970895,
    0.02065439019027026680794212,
    0.02111585199536373124165436,
    0.02158368653362269615798186
  };
  if (!evaluate_and_compare(
        [](int i) { return i * 0.001; },
        [](double x) { return EH::__lagrangian_special2(x); },
        result1,
        sizeof(result1) / sizeof(result1[0]),
        8e-16)) {
    return false;
  }

  static const double result2[] = {
    0.1008275209118719944244002,
    0.09373883714509216768682918,
    0.08660393060662854081119753,
    0.07942512230712844690888330,
    0.07220465097710699957095049,
    0.06494467596502181507863357,
    0.05764728003739881530859239,
    0.05031447208341101704286656,
    0.04294818972638048075842287,
    0.03555030184472078264051480,
    0.02812261100486636931908388,
    0.02066685580874849535600487,
    0.01318471315837736554051857,
    0.005677800440078576026562793,
    -0.001852322369089298832403815,
    -0.009404150670238572835761753,
    -0.01697623328441667623010633,
    -0.02456717056481165902330234,
    -0.03217561257702655727298340,
    -0.03980025734510066079186937,
    -0.04743984916100916590432464,
    -0.05509317695542845669024117,
    -0.06275907272761174155210628,
    -0.07043641003227846661158986,
    -0.07812410252148039051874432,
    -0.08582110253946705516234588,
    -0.09352639976863329596385581,
    -0.1012390199246911239900569,
    -0.1089580234992675420133451,
    -0.1166825045481884267289107,
    -0.1244115895237663508359791,
    -0.1321444361494669912576190,
    -0.1398802323353844579689904,
    -0.1476181951330103880811504,
    -0.1553575697278349073807156,
    -0.1630976284683695084506222,
    -0.1708376699302324872343313,
    -0.1785770180139867874323324,
    -0.1863150210754679042637801,
    -0.1940510510873858850771154,
    -0.2017845028310304312879768,
    -0.2095147931169516582842486,
    -0.2172413600335312172556326,
    -0.2249636622223992403330033,
    -0.2326811781796919570893346,
    -0.2403934055821828689661343,
    -0.2480998606373570840227387,
    -0.2558000774565338353703628,
    -0.2634936074501763624252567,
    -0.2711800187445612558313339,
    -0.2788588956190110868014753,
    -0.2865298379629246927423436,
    -0.2941924607518689069211884,
    -0.3018463935420238344620299,
    -0.3094912799823010240771955,
    -0.3171267773434800985294849,
    -0.3247525560637346205425616,
    -0.3323682993099422180527563,
    -0.3399737025541973061986953,
    -0.3475684731649671556178569,
    -0.3551523300123535991988994,
    -0.3627250030869433735002709,
    -0.3702862331317499869544798,
    -0.3778357712867691243522503,
    -0.3853733787456879647858931,
    -0.3928988264243064362759537,
    -0.4004118946402453819502739,
    -0.4079123728035328963140501,
    -0.4154000591176757314476036,
    -0.4228747602908376966811105,
    -0.4303362912567614053944734,
    -0.4377844749050835822419297,
    -0.4452191418207074556755925,
    -0.4526401300319085457278439,
    -0.4600472847668624364289928,
    -0.4674404582182949160511145,
    -0.4748195093159661959134900,
    -0.4821843035067117983742819,
    -0.4895347125417731547832746,
    -0.4968706142711609918186903,
    -0.5041918924448042263547497,
    -0.5114984365202463507434961,
    -0.5187901414766601874624305,
    -0.5260669076349604391968098,
    -0.5333286404838016717301521,
    -0.5405752505112572560862699,
    -0.5478066530419823762377433,
    -0.5550227680796714918839800,
    -0.5622235201546276443185464,
    -0.5694088381762677187771818,
    -0.5765786552903942399448410,
    -0.5837329087410704891129134,
    -0.5908715397369417019948225,
    -0.5979944933218508451921432,
    -0.6051017182496029861171976,
    -0.6121931668627375748008031,
    -0.6192687949751730550555104,
    -0.6263285617585931251811294,
    -0.6333724296324486827062765,
    -0.6404003641574540211463411,
    -0.6474123339324602067006294,
    -0.6544083104945927561849366,
    -0.6613882682225447709875538,
    -0.6683521842429205618574022,
    -0.6753000383395285320281509,
    -0.6822318128655256774431884,
    -0.6891474926583195183209100,
    -0.6960470649571366014021215,
    -0.7029305193231699121418021,
    -0.7097978475622206158218647,
    -0.7166490436497525108402262,
    -0.7234841036582804308478373,
    -0.7303030256870165793440340,
    -0.7371058097937014250056111,
    -0.7438924579285483324471642,
    -0.7506629738702335551544043,
    -0.7574173631638655787043200,
    -0.7641556330608700766400341,
    -0.7708777924607289319119731,
    -0.7775838518545138868988789,
    -0.7842738232701574178168252,
    -0.7909477202194053878182083,
    -0.7976055576463979201593919,
    -0.8042473518778267512414836,
    -0.8108731205746190757571796,
    -0.8174828826850995851534704,
    -0.8240766583995840285876647,
    -0.8306544691063591948561035,
    -0.8372163373490057266597658,
    -0.8437622867850216371965105,
    -0.8502923421457058055068077,
    -0.8568065291972620832359263,
    -0.8633048747030859534182881,
    -0.8697874063871969433701562,
    -0.8762541528987812105527679,
    -0.8827051437778098940269844,
    -0.8891404094216999564827063,
    -0.8955599810529853343473997,
    -0.9019638906879672676519506,
    -0.9083521711063136985933356,
    -0.9147248558215786094600778,
    -0.9210819790526131181015350,
    -0.9274235756958410636969282,
    -0.9337496812983726984358426,
    -0.9400603320319309530319150,
    -0.9463555646675655668827972,
    -0.9526354165511311682453779,
    -0.9588999255795061570565091,
    -0.9651491301775299839964292,
    -0.9713830692756371350261247,
    -0.9776017822881668218551505,
    -0.9838053090923280465013013,
    -0.9899936900078003531420223,
    -0.9961669657769512036507083,
    -1.002325177545651515348583,
    -1.008468366844671481343943,
    -1.014596575571639356105316,
    -1.020709845973546432325848,
    -1.026808220629781960358365,
    -1.032891742435682269184006,
    -1.038960454586578838646987,
    -1.045014400562330547145613,
    -1.051053624112325777693557,
    -1.057078169240940508812750,
    -1.063088080193438945625771,
    -1.069083401442303661297211,
    -1.075064177673982620126636,
    -1.081030453776040841598279,
    -1.086982274824704840004717,
    -1.092919686072788337326283,
    -1.098842732937988098291544,
    -1.104751460991539076377400,
    -1.110645915947218388325866,
    -1.116526143650687952939217,
    -1.122392190069165937832970,
    -1.128244101281417455830437,
    -1.134081923468055241113926,
    -1.139905702902141314434079,
    -1.145715485940080916936428,
    -1.151511319012800253797594,
    -1.157293248617199842166172,
    -1.163061321307875503158130,
    -1.168815583689099275136163,
    -1.174556082407052755469696,
    -1.180282864142305600680541,
    -1.185995975602532130572755,
    -1.191695463515459190859696,
    -1.197381374622038631163529,
    -1.203053755669837951291733,
    -1.208712653406642858602463,
    -1.214358114574265663259869,
    -1.219990185902553616447972,
    -1.225608914103591469347127,
    -1.231214345866092698063228,
    -1.236806527849974001913138,
    -1.242385506681107839680508,
    -1.247951328946247920828527,
    -1.253504041188122716348698,
    -1.259043689900692197090447,
    -1.264570321524563146202978,
    -1.270083982442558526870734
  };
  if (!evaluate_and_compare(
        [](int i) { return (i + 100) * 0.01; },
        [](double x) { return EH::__lagrangian_special2(x); },
        result2,
        sizeof(result2) / sizeof(result2[0]),
        8.1e-14)) {
    return false;
  }

  static const double result3[] = {
    0,
    0.02158368653362269615798186,
    0.1003440792063249937678914,
    0.2033614252238973770626031,
    0.2771242272362825700402149,
    0.3091540929018506963460622,
    0.3057034633138983277406588,
    0.2761390620828977878490501,
    0.2284497419050731085263673,
    0.1685963368646018698446581,
    0.1008275209118719944244002,
    0.02812261100486636931908388,
    -0.04743984916100916590432464,
    -0.1244115895237663508359791,
    -0.2017845028310304312879768,
    -0.2788588956190110868014753,
    -0.3551523300123535991988994,
    -0.4303362912567614053944734,
    -0.5041918924448042263547497,
    -0.5765786552903942399448410,
    -0.6474123339324602067006294,
    -0.7166490436497525108402262,
    -0.7842738232701574178168252,
    -0.8502923421457058055068077,
    -0.9147248558215786094600778,
    -0.9776017822881668218551505,
    -1.038960454586578838646987,
    -1.098842732937988098291544,
    -1.157293248617199842166172,
    -1.214358114574265663259869,
    -1.270083982442558526870734,
    -1.324517357557279647254850,
    -1.377704106703403539874176,
    -1.429689110107222135969126,
    -1.480516021481693757827694,
    -1.530227108994751759567637,
    -1.578863156743041005922426,
    -1.626463411315737827418526,
    -1.673065561779089063365445,
    -1.718705744230510826273495,
    -1.763418564200208225808830,
    -1.807237131792880133801770,
    -1.850193105690594507137097,
    -1.892316743075394446151873,
    -1.933636953247322041196428,
    -1.974181353263189575510916,
    -2.013976324343303577879414,
    -2.053047068117433081758504,
    -2.091417662030330212093894,
    -2.129111113418349288384099,
    -2.166149411915376930843272,
    -2.202553579958503046577535,
    -2.238343721249414930536377,
    -2.273539067092387003548219,
    -2.308158020578608835883043,
    -2.342218198623042148753376,
    -2.375736471886818345546583,
    -2.408729002637556777539783,
    -2.441211280613595527667764,
    -2.473198156967303723423872,
    -2.504703876368420255828303,
    -2.535742107351545363562605,
    -2.566325970993129418176629,
    -2.596468068003049858572779,
    -2.626180504314526734877357,
    -2.655474915253998920383084,
    -2.684362488369899166550081,
    -2.712853984996206763958279,
    -2.740959760623360318825642,
    -2.768689784145686297200249,
    -2.796053656051022257773567,
    -2.823060625614747789673985,
    -2.849719607157026086493957,
    -2.876039195418737523355088,
    -2.902027680108376669583867,
    -2.927693059669101518627310,
    -2.953043054312178234057667,
    -2.978085118360261847045684,
    -3.002826451941295167919319,
    -3.027274012071294215437777,
    -3.051434523161916351651908,
    -3.075314486986473287863218,
    -3.098920192135950471921637,
    -3.122257722994621690289056,
    -3.145332968262997212776816,
    -3.168151629054109435840737,
    -3.190719226587515616946673,
    -3.213041109503876834438649,
    -3.235122460821549753458787,
    -3.256968304555297288916538,
    -3.278583512015980194435114,
    -3.299972807808928566157258,
    -3.321140775547605076137855,
    -3.342091863298155544848782,
    -3.362830388769492592021144,
    -3.383360544262670197505825,
    -3.403686401392476944334738,
    -3.423811915593399640419846,
    -3.443740930421383300634788,
    -3.463477181662134726223557,
    -3.483024301256081973119441
  };
  if (!evaluate_and_compare(
        [](int i) { return i * 0.1; },
        [](double x) { return EH::__lagrangian_special2(x); },
        result3,
        sizeof(result3) / sizeof(result3[0]),
        5e-15)) {
    return false;
  }

  return true;
}

bool test_euler_heisenberg(void) {
  /* Manually testing partial results. */
  std::pair<double, double> result = EH::__lagrangian_ab_helper(0.14, 0.063);
  if (!compare_eq_rel(result.first, 2.76727232e-6, 1e-5, 0.))
    return false;
  if (!compare_eq_rel(result.second, 3.094011654258106105473541e-25, 4e-4, 0.))
    return false;

  /* Compare with the lowest order Lagrangian. */
  for (int k = 0; k < 100; ++k) {
    double F = random_double(0.0, 1e-6 * sqr(PHY_Bc));
    double G = random_double(0.0, 1e-6 * sqr(PHY_Bc));
    if (!compare_eq_rel(EH::lagrangian_real(F, G),
                        QED::lagrangian_lowest_order(F, G),
                        1e-5, 0.)) {
      return false;
    }
  }

  /* Compare with the magnetic-only general Lagrangian. */
  for (int k = 0; k < 100; ++k) {
    double F = random_double(0.0, 1e-6 * sqr(PHY_Bc));
    double G = 0.0;
    if (!compare_eq_rel(EH::lagrangian_real(F, G),
                        magnetic_only_lagrangian(F, G),
                        1e-7, 0.)) {
      fprintf(stderr, "k=%d F=%lg G=%lg\n",
          k, F / sqr(PHY_Bc), G / sqr(PHY_Bc));
      return false;
    }
  }

  /* Compare with the lowest order lambdas. */
  for (int k = 0; k < 100; ++k) {
    // Low-limit.
    double F = random_double(0.0, 1e-10 * sqr(PHY_Bc));
    double G = random_double(0.0, 1e-10 * sqr(PHY_Bc));
    std::pair<double, double> lambdas = qed_metric_correction_lambda(
        [](auto F, auto G) { return EH::lagrangian_real(F, G); }, F, G);
    bool ok1 = compare_eq_rel(lambdas.first, QED::lambda1, 1e-5, 0.);
    bool ok2 = compare_eq_rel(lambdas.second, QED::lambda2, 1.2e-3, 0.);
    if (!ok1 || !ok2) return false;
  }

  return true;
}
#endif
