<map version="0.8.1">
<!-- To view this file, download free mind mapping software FreeMind from http://freemind.sourceforge.net -->
<node CREATED="1269766250289" ID="Freemind_Link_531478983" MODIFIED="1269779361504" TEXT="Density">
<node CREATED="1269766262263" ID="_" MODIFIED="1269780587454" POSITION="right" TEXT="read 1)&#xa; MT3DRhoFlag, &#xa; MFNADVFD,&#xa; NSWTCPL,&#xa; IWTABLE">
<node CREATED="1269779716447" ID="Freemind_Link_678409351" MODIFIED="1269779765112" TEXT="if NSWTCPL&lt;&gt;0">
<node CREATED="1269779570835" ID="Freemind_Link_485768584" MODIFIED="1269779581607" TEXT="read 3) DensCrit"/>
</node>
<node CREATED="1269766281696" ID="Freemind_Link_611804531" MODIFIED="1269780912494" TEXT="if MT3DRhoFlag&gt;=0, read 4)">
<node CREATED="1269779970281" ID="Freemind_Link_35696591" MODIFIED="1269781019020" TEXT="if MT3DRhoFlag &gt;0)&#xa;MT3DRhoFlag = MT3D species for density computation ">
<node CREATED="1269780969593" ID="Freemind_Link_1389718736" MODIFIED="1269781002814" TEXT="Read 2)&#xa;DenseRef, dRhodC"/>
</node>
<node CREATED="1269769651986" ID="Freemind_Link_477761451" MODIFIED="1269780567766" TEXT="if MT3DRhoFlag=0">
<node COLOR="#ff0000" CREATED="1269769609819" ID="Freemind_Link_1260957548" MODIFIED="1269770509696" TEXT="FOR EACH STRESS PERIOD">
<node CREATED="1269767807227" ID="Freemind_Link_1631965123" MODIFIED="1269780644068" TEXT="read 6)&#xa;InDense">
<node CREATED="1269781960591" ID="Freemind_Link_1450874475" MODIFIED="1269781969587" TEXT="if InDense &lt; 0">
<node CREATED="1269781970503" ID="Freemind_Link_1032283353" MODIFIED="1269782012880" TEXT="Reuse from previous stress period&#xa;if first stress period set all cells to DenseRef"/>
</node>
<node CREATED="1269782032220" ID="Freemind_Link_1079075432" MODIFIED="1269782039375" TEXT="if InDense = 0">
<node CREATED="1269782040612" ID="Freemind_Link_1496951574" MODIFIED="1269782052672" TEXT="Set all cells to DenseRef"/>
</node>
<node CREATED="1269767827740" ID="Freemind_Link_861165419" MODIFIED="1269780656704" TEXT="if InDense &gt;0 read 7)">
<node CREATED="1269767834772" ID="Freemind_Link_269004612" MODIFIED="1269780675265" TEXT="Dense(NROW,NCOL,Layer1)"/>
<node CREATED="1269767862964" ID="Freemind_Link_813392645" MODIFIED="1269780682849" TEXT="Dense(NROW,NCOL,Layer2)"/>
<node CREATED="1269767876301" ID="Freemind_Link_1727037609" MODIFIED="1269767876301" TEXT=""/>
<node CREATED="1269767878053" ID="Freemind_Link_789615689" MODIFIED="1269780691033" TEXT="Dense(NROW,NCOL,NLAY)"/>
</node>
<node CREATED="1269782087133" ID="Freemind_Link_1465239872" MODIFIED="1269782094241" TEXT="if InDense = 2">
<node CREATED="1269782094613" ID="Freemind_Link_1924980396" MODIFIED="1269782131570" TEXT="Values read are concentrations&#xa;the will be converted to densities&#xa;by Seawat"/>
</node>
</node>
</node>
</node>
</node>
<node CREATED="1269766289464" ID="Freemind_Link_852213357" MODIFIED="1269780509149" TEXT="if MT3DRhoFlag = -1">
<node CREATED="1269780128274" ID="Freemind_Link_937382233" MODIFIED="1269780194193" TEXT="read 4a)&#xa;DenseRef, dRhodPrHd, PrHdRef"/>
<node CREATED="1269766368090" ID="Freemind_Link_1379144046" MODIFIED="1269780411524" TEXT="read 4b)&#xa;NSRhoEOS">
<node CREATED="1269766444973" ID="Freemind_Link_81369224" MODIFIED="1269780427530" TEXT="if NSRhoEOS&gt;0 read 4c)&#xa;NSRhoEOS = # of species to&#xa;include in the&#xa;density computation">
<node CREATED="1269766476138" ID="Freemind_Link_1759568721" MODIFIED="1269781229188" TEXT="MTRhoSpec(1), DRhoDc(1), cRhoRef(1)"/>
<node CREATED="1269766524010" ID="Freemind_Link_1245394456" MODIFIED="1269781238956" TEXT="MTRhoSpec(2), DRhoDc(2), cRhoRef(2)"/>
<node CREATED="1269766527594" ID="Freemind_Link_366743975" MODIFIED="1269781246972" TEXT="MTRhoSpec(...), DRhodc(...), cRhoRef(...)"/>
<node CREATED="1269766534362" ID="Freemind_Link_1571136472" MODIFIED="1269781258573" TEXT="MTRhoSpec(N), DRhodc(N), cRhoRef(N)"/>
</node>
</node>
</node>
</node>
<node CREATED="1269767517179" ID="Freemind_Link_1018212425" MODIFIED="1269779507058" POSITION="right" TEXT="read 2)&#xa;DensMin, DensMax"/>
<node CREATED="1269779585795" ID="Freemind_Link_356473620" MODIFIED="1269779595831" POSITION="right" TEXT="read 5) FirstDt"/>
</node>
</map>
