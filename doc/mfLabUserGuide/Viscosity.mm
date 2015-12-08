<map version="0.8.1">
<!-- To view this file, download free mind mapping software FreeMind from http://freemind.sourceforge.net -->
<node CREATED="1269766250289" ID="Freemind_Link_531478983" MODIFIED="1269766260707" TEXT="Viscosity">
<node CREATED="1269766262263" ID="_" MODIFIED="1269781370756" POSITION="right" TEXT="read 1)&#xa;MT3DMuFlag">
<node CREATED="1269766281696" ID="Freemind_Link_611804531" MODIFIED="1269781425692" TEXT="if MT3DMuFlag&gt;=0, read 3)">
<node CREATED="1269766299920" ID="Freemind_Link_1238779524" MODIFIED="1269781593319" TEXT="if MT3DMuFlag&gt;0, read 3)&#xa;MT3DMuFlag = MT3D species for Temp">
<node CREATED="1269781483286" ID="Freemind_Link_1192245695" MODIFIED="1269781494987" TEXT="ViscRef, dMudC, cMuRef"/>
</node>
<node CREATED="1269769651986" ID="Freemind_Link_477761451" MODIFIED="1269781614414" TEXT="if MT3DMuFlag=0">
<node COLOR="#ff0000" CREATED="1269769609819" ID="Freemind_Link_1260957548" MODIFIED="1269770509696" TEXT="FOR EACH STRESS PERIOD">
<node CREATED="1269767807227" ID="Freemind_Link_1631965123" MODIFIED="1269768628312" TEXT="read 4)&#xa;InVisc">
<node CREATED="1269782169270" ID="Freemind_Link_601123935" MODIFIED="1269782180953" TEXT="if InVisc &lt; 0">
<node CREATED="1269782297271" ID="Freemind_Link_348163208" MODIFIED="1269782325196" TEXT="Values of previous period are resued&#xa;if first stress period, then values are set to ViscRef"/>
</node>
<node CREATED="1269782183294" ID="Freemind_Link_1724617683" MODIFIED="1269782187849" TEXT="if InVisc = 0">
<node CREATED="1269782282880" ID="Freemind_Link_1006750763" MODIFIED="1269782294899" TEXT="Values are set to ViscRef"/>
</node>
<node CREATED="1269767827740" ID="Freemind_Link_861165419" MODIFIED="1269769186701" TEXT="if InVisc &gt;0 read 5)">
<node CREATED="1269767834772" ID="Freemind_Link_269004612" MODIFIED="1269767862336" TEXT="Visc(NROW,NCOL,Layer1)"/>
<node CREATED="1269767862964" ID="Freemind_Link_813392645" MODIFIED="1269767875049" TEXT="Visc(NROW,NCOL,Layer2)"/>
<node CREATED="1269767876301" ID="Freemind_Link_1727037609" MODIFIED="1269767876301" TEXT=""/>
<node CREATED="1269767878053" ID="Freemind_Link_789615689" MODIFIED="1269767891185" TEXT="Visc(NROW,NCOL,NLAY)"/>
</node>
<node CREATED="1269782191781" ID="Freemind_Link_1596748363" MODIFIED="1269782200841" TEXT="if InVisc = 2">
<node CREATED="1269782210446" ID="Freemind_Link_523367828" MODIFIED="1269782240376" TEXT="Values read are concentrations&#xa;Seawat will convert them to viscosities"/>
</node>
</node>
</node>
</node>
</node>
<node CREATED="1269766289464" ID="Freemind_Link_852213357" MODIFIED="1269781629286" TEXT="if MT3DMuFlag = -1">
<node CREATED="1269766368090" ID="Freemind_Link_1379144046" MODIFIED="1269781663421" TEXT="read 3b)&#xa;NSMuEOS, MuTempOpt">
<node CREATED="1269766444973" ID="Freemind_Link_81369224" MODIFIED="1269781682656" TEXT="if NSMuEOS&gt;0 read 3c)&#xa;NSMuEOS = # of species to&#xa;include in the&#xa;viscosity computation">
<node CREATED="1269766476138" ID="Freemind_Link_1759568721" MODIFIED="1269781698864" TEXT="MTMuSpec(1), dMudC(1), cMuRef(1)"/>
<node CREATED="1269766524010" ID="Freemind_Link_1245394456" MODIFIED="1269781715920" TEXT="MTMuSpec(2), dMudC(2), cMuRef(2)"/>
<node CREATED="1269766527594" ID="Freemind_Link_366743975" MODIFIED="1269781725353" TEXT="MTMuSpec(...), dMudC(...), cMuRef(...)"/>
<node CREATED="1269766534362" ID="Freemind_Link_1571136472" MODIFIED="1269781737393" TEXT="MTMuSpec(N), dMudC(N), cMuRef(N)"/>
</node>
<node CREATED="1269766538698" ID="Freemind_Link_1082204747" MODIFIED="1269781839832" TEXT="if MuTempOpt &gt;0, read 3d)&#xa;MuTempOpt = vicosity eq #">
<node CREATED="1269766631255" ID="Freemind_Link_1579881711" MODIFIED="1269781761916" TEXT="if MuTempOpt =3&#xa;Visosity eq 20">
<node CREATED="1269767100566" ID="Freemind_Link_1416186615" MODIFIED="1269781789706" TEXT="read MTMuTempSpec +2 coefficients"/>
</node>
<node CREATED="1269766623599" ID="Freemind_Link_899973927" MODIFIED="1269768997108" TEXT="if MUEMPTOPT=2&#xa;viscosity eq 19">
<node CREATED="1269767083670" ID="Freemind_Link_172348551" MODIFIED="1269781798530" TEXT="read MTMuTempSpec +5 coefficients"/>
</node>
<node CREATED="1269766612150" ID="Freemind_Link_791524864" MODIFIED="1269769024417" TEXT="if MUEMPTOP=1&#xa;viscosity eq 18">
<node CREATED="1269767053278" ID="Freemind_Link_921425731" MODIFIED="1269781808796" TEXT="read MTMuTempSpec +4 coefficients"/>
</node>
</node>
<node CREATED="1269766707605" ID="Freemind_Link_271974140" MODIFIED="1269781860219" TEXT="if (NSMuEOS=0 &amp; MuTempOpt=0) ">
<node CREATED="1269766729269" ID="Freemind_Link_1256028593" MODIFIED="1269766748305" TEXT="viscosity=ViscRef"/>
</node>
</node>
<node CREATED="1269767570749" ID="Freemind_Link_181856240" MODIFIED="1269768693155" TEXT="read 3a)&#xa;ViscRef"/>
</node>
</node>
<node CREATED="1269767517179" ID="Freemind_Link_1018212425" MODIFIED="1269768612192" POSITION="right" TEXT="read 2)&#xa;ViscMin, ViscMax"/>
</node>
</map>
