__author__ = 'ilap'

from SequenceCompareLib import *

def getMaxOfNodes (graph, row, col):

    row_len = len(graph)-1
    col_len = len (graph[0,:])-1

    maxv = -sys.maxint

    to_row = row
    to_col = col

    if col_len > col:
        tv = graph[row,col+1]
        if maxv < tv:
            maxv = tv
            to_col = col + 1


    if row_len > row:
        tv = graph[row+1,col]
        if maxv < tv:
            maxv = tv
            to_row = row + 1

    if row_len > row and col_len > col:
        tv = graph[row+1,col+1]
        if maxv < tv:
            maxv = tv
            to_row = row + 1
            to_col = col + 1

    return [(row, col), (to_row, to_col)]
    #if row_len


def getMiddleEdges (alligned_graph):

    col_idx = int(len (alligned_graph[0,:])/2-1)

    tcol = alligned_graph[:,col_idx]
    max_value = tcol.max()
    rows = [i for i, j in enumerate(tcol) if j == max_value]
    print rows

    for row_idx in rows:
        print getMaxOfNodes (alligned_graph, row_idx, col_idx)

def middleNodeGraph (graph, top, bottom, left, right):
    middle_idx = int ((right-left-1)/2)
    tcol = b[top:bottom,middle_idx]
    max_value = tcol.max()

    return list (tcol).index(max_value)

    #rows = [i for i, j in enumerate(tcol) if j == max_value]

def middleNodeGraph (graph, top, bottom, left, right):
    middle_idx = int ((right-left-1)/2)
    tcol = b[top:bottom,middle_idx]
    max_value = tcol.max()

    return list (tcol).index(max_value)


def middleEdge (graph, top, bottom, left, right):

    col_idx = int(len (alligned_graph[0,:])/2-1)

    tcol = alligned_graph[:,col_idx]
    max_value = tcol.max()
    rows = [i for i, j in enumerate(tcol) if j == max_value]
    print rows

    for row_idx in rows:
        print getMaxOfNodes (alligned_graph, row_idx, col_idx)


def linearSpaceAlignment (v, w, top, bottom, left, right):
    if left == right:
        return "D"
    if top == bottom:
        return "C"

    middle = int ((left+right)/2)

    middle_node = middleNode (top, bottom, left, right)
    middle_edge = middleEdge (top, bottom, left, right)

    linearSpaceAlignment(top, middle_node, left, middle)
    print middle_edge

    if middle_edge == "r" or middle_edge == "m":
        middle +=1
    if middle_edge == "d" or middle_edge == "m":
        middle_node += 1

    linearSpaceAlignment(v, w, middle_node, bottom, middle, right)



r="PLEASANTLY"
c="MEASNLY"

rh="WTGWGIFKEGIPCVHIMREHQHLQLAWGIESMIAYAWHDQCIPWCVRYMAAGNGLYDNYKNICRGTKTPITVENNYGPENMVKRPLINVPEAGNMNFTFMWIVYWTMWCICKIRPMPITFNCYMLAEQLNEWVVGGACNGMHQVSVKSVYIVHPPLKHGCERTCYKRFSWAGLQAEKEMVWAQMAVSYGIFHRHGWHQWNDLTFSLRCTVVQYWYAFPRMCHETVHALEMNQYLQEIGYLQHEERYGAYMLWMLFECEATEFSYHTRNCMVDSHFYELCMMSHCFEVQFMYQHKFYQRGDFQFGGWPLIEYGIAKCSRCRPPTSMWYYSTIFPWIDPIYAKAVSSPSSCWTYATYWKRTLYPQPNIWEMFHMGMESMYLWLVTCNKMKTKANVAPFQNNFLPVGHDTPWNGKSFGKNFPNKWDADEKACLWYVTKMMQVYGTAVDLVEMSNYMFSYSRSMPRYETTGYVKPWKWHHMMIEAMNIYYVDELLRKLWYFYVVYETHYPPFSIIYDFSNNEWKAEIIWKVYGINSKDWKAVWWEQWNTYWAQKPKGMWRLCVPVQGFTGKTQYCKFNEVFHFHSIRPCDLREPDNVTHPINREFIYWGSFQPPPAMKHAPVNSCFYMFRCCSSKVHGTRGWCIYSENMCADCLTEKLVDIVIWCIDTCYPAACTHYSEIFCAVVHSSTPTSWLQSWRLSYLWHRWFGVTGVDHNGLWNDRIHDPNGQIKTHDYAQHSYKNYQYALQLLPFKWGRDNFMTLIHWFIKCDQEVMVEMHHKTFRKMVMKFEKRPMYTQKRASQYIRKFPKWTNPQFTSWQDHNCFMAHLSLRHNHGPKDDRINWKNMDAPYVKCCGTVSVSAILLMWPNAKEMLQFHLWSRCNYMRDGLWGIWGQKPMVGPWIEMYNKDCPIVGNWQYKANNIIVYVLCKVDKGQETNESSEKIWSMLFSYPLLCTKRAQKMNTCEWAMSESKRQHLDFCKDSMTTMCNCYYWTKGQEFTFMVYYRYYSCHEKESQMPLRRYEWVENSTHDWHGYGPLNCK"
ch="GYQNAEVSHNDFCEDLKTVRHHNRMLGFVTAHHYSYPEKPCEFMTYTYRMSNYLVWVTLDFDNQGKNPHHRKDGKDFWDFGDWNWRYCDCCYVDAGFFKTRMCCCNPYHDWFHFCACPFAPVKTTQPIWKFWISMCNKVVYHSKNIIIDRILGEAGRPAKVKFELTRGQAKLFQLQCIFAHVHTHHGLWGNMSVVGEVHIRTTWQFACWIIFYWRTNPFMNWGAYPSQWDQATIEQDTKRRSCRTKSELNRYCDYCHQHKVSKEGCWYIHQPCTDHAHMIKWEDQQPGWFATFECGAHDIRSIYKSQPMGQKMERRDAYPIWQQDIKYHNQRCDMYVWSEMKEERSAWACHMWSVYPIGDATVDWHQMYLTLVRSEEPLQWNNMFRLMVKYHWICELHFNCGFIQTVDRYEKNSGLRNSKIAKYLHTKSKRSAMHWVIDEPPNSTPPRFIWYMYNSRSMPRYETTGYVKPWKWHDKRNVGHMGLRIEAMNIYYVIELLRLLKEFYVVYHMAHDYVPFSIIYFFLNPEWKAEIYWKVYGINWVGRCCSSVDWKAVWWEQWNTYWAQKPKGMWCLVPVLGFTGKTQYCKFNEVFHFCSIRPCWHSNHVWESYQRALRDKHHCNYDLGLVASRVRPCEDNQYTKCQNENADLFCFKFMRFIDATERSSKDFHFVAEIASQIIVQFKVVMRWDGWVNKPFAGMMHFCDAQHPDEYRRWSYKIDKKYMDRGRNDYADYQGRFHQGEVDPVWSYDNGERNCRNIWPRMNKYMFASEHLSMNMNPLDSGKRYQLEKPDNFRLADRVTWRTPTMYGGMERNFEGCAQLMSADLQEWDKARITCPMLHPNPDAYSPYRYQQIPDRVYLFEKIHDHWGEGVSKRKIRHGKRSWGFMLLHNQTNLAIQMENCGDCFVIMPNLLTNGLFLPVEFDLQTDMPAMHDHNHRVIPLNDGTWEQNIHLICVQHYYSAVPQMTNPVTKSFRKWAGHMICAFKSPHANAQIYFYAKGHNMCLIHWCPLARKPAMNMYLFCKYSRVNVIMMWGMFPEVYKVDELTFLGHL"

ro = 5
(b, row, col, matrix) =  globalAlignmentGraph(r, c, ro)
print b

top=0
left=0
bottom = len (r)+1
right = len (c)+1

print middleNode(b, top, bottom, left, right)

print getMiddleEdges()

#print rows



#print getMiddleEdges(b)

#row_idx = list(middle_col).index()
#print row_idx
