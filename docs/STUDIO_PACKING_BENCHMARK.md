# Kaynakları eşlenmiş MXene–polimer paketleme karşılaştırması

**Sistem oluştur** ekranında dört paketleyici kullanılabilir. `builtin` önceki PolyChemMD yöntemidir. `slit` uzun bir rijit zinciri iki MXene tabakasının arasına yerleştirmek için geliştirilmiştir. Diğer iki seçim [MolCrafts MolPack](https://docs.molcrafts.org/molpack/) ve [Packmol](https://doi.org/10.1002/jcc.21224) motorlarını çalıştırır. Dört motor da aynı atom kimliği ve koordinat denetiminden geçer. `packing_seed` yerleştirme için kullanılır. Bu alan boş bırakılırsa önceki tariflerin `seed` değeri kullanılır.

## Yerleşim ve denetim

Yeni yöntemde zincirin ana eksenleri merkezlenmiş atom koordinatlarının tekil değer ayrışımından bulunur. Her iki denemeden birinde en ince eksen tabaka normaline döndürülür. Düzlem içi azimut ve iki yüz yönünden biri rastgele seçilir. Diğer denemede izotropik bir üç boyutlu rotasyon kullanılır. Seçilen yönelimde bütün atomları kutuya sığdıracak ötelenme aralığı her eksende doğrudan hesaplanır. Boş olmayan aralıktan uniform öteleme önerilir. Zincirin iç bağları ve kiralitesi korunur. İyonlar ve sular rijit yerleştirilir. Bir yerleşim çıkmaza girerse en fazla iki yeni başlangıç yapılır. Rapor başarısız başlangıçları ve toplam yerleştirme denemesi sayısını saklar.

Bu yöntem başlangıç konformasyonunun tabakaya göre yönelimini etkileyebilir. Adsorpsiyon, su düzeni veya zincir yönelimi hakkında çıkarım için bağımsız başlangıçlardan dengeye yaklaşmanın gözlenebilir düzeyinde sınanması gerekir. Katmanlar arasında sığmayan bir izole zincirin iç torsiyonları paketlemede değiştirilmez.

Bütün motorlar için ilk fiziksel olabilirlik testi, şablonun en büyük atomlar arası uzaklığını erişilebilir kutunun köşegeniyle karşılaştırır. Ek olarak, merkezlenmiş koordinat kovaryansının en küçük özdeğeri λ ise, her doğrultuda şablon genişliği en az 2√λ olur. Bu iki eşitsizliği geçmek çözümü garanti etmez. Çıktı denetimi tüm molekül içi mesafeleri, dört veya daha fazla atomlu şablonlarda yansımasız rijit dönüşümü, sabit yüzeyi, kutu sınırını, XY periyodik görüntüleri ve farklı moleküller arasındaki en kısa teması inceler.

## Eşlenmiş karşılaştırma

| Örnek | Polimer ve sıvı | Hücre |
| --- | --- | --- |
| `aqueous` | PEO diol, DP8, 1 zincir, 0,997 g/cm³ nominal su, tuzsuz | Tek tabaka. Su yüksekliği 22 Å, sıvının üstünde 12 Å boşluk; üst sınır açık yüzeyin hemen yanında değildir. |
| `confined` | %75 hidrolizli PVA/VAc, DP12, 1 zincir, kümeli asetatlar, 0,85 g/cm³ nominal su, tuzsuz | İki tabaka. İç uç atom çekirdekleri arasındaki başlangıç açıklığı 9 Å. |

Her örnekte aynı 6 × 4 bazal üst hücre ve aynı tohumla üretilen tek bir zincir/yüzey geometrisi bütün motorlara verilir. Her yerleştirme tohumuyla dört motor ayrı başlatılır. Gerçekleşen su sayısı ve atom sayısı her başarılı koşuda kaydedilir. Ön deneme sürecinde, tek tabakalı kutuda üst hücre sınırını sıvının üst yüzüne yerleştirmek bazı LAMMPS enerji minimizasyonlarında yüzeyden atom kaybı verdi. Sabit sıvı hacmine eklenen boşluk bu nedenle tarifte açıkça yer alır.

Geometrik kabul eşiği 1,80 Å, ortak üretim toleransı bunun 0,035 Å üstü ve kutu sınır payı 0,15 Å'dır. XYZ okuma hassasiyeti için kabul toleransı 0,002 Å'dır. Tüm araçlara tohum başına 60 s alt süreç bütçesi ve sayısal motorlara bir iş parçacığı sınırı verilir. Yerel yöntem bir kopya için 3.500 merkez/yönelim denemesi, harici motorlar 250 döngü üst sınırı kullanır. Bu değerler farklı algoritmaların kendi çaba birimleridir. Yerel `slit` motorunun yeniden başlatmaları aynı zaman bütçesinden düşülür. Zaman ölçümü paketleme ve bağımsız son kontrolleri kapsar. Python başlangıcını da kapsayan süre ayrı sütundadır. Tohumlar arasında motor sırası döndürülür.

`evaluate_forces=true` seçildiğinde, başarılı her geometri için aynı bileşimde ve aynı kuvvet alanı kurallarında yeni topoloji oluşturulur. LAMMPS sabit substratlı kutuda 8 Å kesme ile en fazla 60 eşlenik gradyan iterasyonu uygular. Mobil atomlardaki en büyük kuvvet normu (kcal mol⁻¹ Å⁻¹), tam hücre potansiyel enerjisi (kcal mol⁻¹), gerçekleşen iterasyon sayısı ve girdi SHA-256 değeri saklanır. Enerji minimizasyonunun ayrı zaman bütçesi 150 s'dir. Başarılı bir geometrinin kuvvet hesabı başarısız olabilir. Bu nedenle iki sonuç ayrı sayılır. Tam kutu enerjisi, özellikle farklı örnekler arasında, bağlanma serbest enerjisi olarak yorumlanamaz.

## Çalıştırma ve kayıt

Stüdyo arayüzünden **Paketleyicileri karşılaştır** formunu açın. Aynı protokole Python ortamından da erişilir.

```sh
python -m polychemmd.studio.packing_benchmark --config examples/studio/packing_benchmark.json --output artifacts/packing-benchmark
```

Kuvvet karşılaştırması için yerel IFF su kaynağı, Moltemplate ve LAMMPS hazırlanmış olmalıdır. `packing_benchmark.json` her koşunun sonucunu, reçeteyi, yazılım sürümlerini, başarı/başarısızlık sayısını ve iki motorun birlikte başarılı olduğu eşlerde medyan zaman farkını içerir. `benchmark_trials.csv` tek satır/deneme biçimindedir. Eksik motor ve süre aşımı ayrı durumlardır. [Depoya eklenen gerçek hesap kayıtları](../evidence/studio/packing_benchmark_v2.json).

## Bu sürümde ölçülen sonuçlar

Aşağıdaki 24 deneme 17, 31 ve 53 yerleştirme tohumlarıyla çalıştırıldı. Sulu model 1.292, galeri modeli 1.169 atom içerir. Başarı geometrik denetimlerin tamamından geçmeyi gerektirir. Kuvvet sütunu, kabul edilen ve 60 iterasyona kadar minimize edilen hücrelerin son mobil maksimum kuvvetinin medyanıdır.

| Koşul | Motor | Geometrik başarı | Başarılı işlerde medyan paketleme süresi (s) | Medyan son Fmax (kcal mol⁻¹ Å⁻¹) |
| --- | --- | ---: | ---: | ---: |
| Sulu PEO | PolyChemMD yerleşik | 3/3 | 0,52 | 27,29 |
| Sulu PEO | PolyChemMD galeri | 3/3 | 0,49 | 60,31 |
| Sulu PEO | MolPack | 3/3 | 0,39 | 47,83 |
| Sulu PEO | Packmol | 3/3 | 0,64 | 17,24 |
| 9 Å galeri, PVA/VAc | PolyChemMD yerleşik | 0/3 | — | — |
| 9 Å galeri, PVA/VAc | PolyChemMD galeri | 1/3 | 1,52 | 23,84 |
| 9 Å galeri, PVA/VAc | MolPack | 0/3 | — | — |
| 9 Å galeri, PVA/VAc | Packmol | 2/3 | 1,77 | 48,09 |

Yerleşik motor galeride üç tohumda da zinciri yerleştirirken bütçesini tüketti. Yeni yöntem bütün tohumlarda zinciri sığdırdı. İki tohumda son su moleküllerinin eklenmesi sırasında yeniden başlatmalar da tükendi. Packmol'ün tohum 17 çıktısında zincirin periyodik görüntüsüyle çakışması bağımsız denetimde yakalandı. MolPack'in üç çıktısı kendi yakınsama ölçütünü karşılamadı. Üç tohum bir zincir ve bir yüzey üzerinde koşulmuştur; polimer aileleri için genelleme ve denge değerlendirmesi ilave bağımsız yapılar gerektirir. [Ham denemeler (CSV)](../evidence/studio/packing_benchmark_v2.csv).

## Yeni yerleştirme tohumlarında tekrar

Aynı 9 Å PVA/VAc hücresi ve aynı motor ayarları, ilk karşılaştırmadan sonra seçilen 101, 211 ve 307 tohumlarında tekrar çalıştırıldı. Yeni galeri motoru 2/3, Packmol 3/3 başarıya ulaştı. Eski yerleşik motor ve MolPack 0/3 kaldı. Beş kabul edilen sistemde sabit yüzey ve aynı minimizasyon protokolü ile kuvvet ölçümü de tamamlandı. İki seri birlikte değerlendirildiğinde altı yerleştirme tohumu için karşılık gelen oranlar 3/6, 5/6, 0/6 ve 0/6'dır. Bu tekrar başlangıç yerleşimini değiştirir. Aynı zinciri ve yüzeyi kullandığı için başka polimerlere veya galeri genişliklerine genellemeyi sınamaz.

[Yeniden çalıştırma reçetesi](../examples/studio/packing_benchmark_replication.json) · [Özet ve tüm denemeler](../evidence/studio/packing_confined_replication.json) · [Ham CSV](../evidence/studio/packing_confined_replication.csv)
