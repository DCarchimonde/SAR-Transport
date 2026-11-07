"""Original vector figures from saved IC50 and completed Ki results. No fitting."""
from pathlib import Path
import json,time
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch,FancyArrowPatch
ROOT=Path(__file__).resolve().parents[1];R=ROOT/'reports/sar_resource_v1';O=ROOT/'manuscripts/sar_resource_v1/figures';O.mkdir(parents=True,exist_ok=True)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'svg.fonttype':'none','pdf.fonttype':42,'axes.spines.top':False,'axes.spines.right':False})
start=time.process_time()
fig,ax=plt.subplots(figsize=(7.1,5.1));ax.set(xlim=(0,10),ylim=(0,8.3));ax.axis('off')
def box(x,y,w,h,title,body,color):
 ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle='round,pad=0.13',edgecolor=color,facecolor=color+'12',linewidth=1.3))
 ax.text(x+.12,y+h-.2,title,va='top',weight='bold',fontsize=10,color=color)
 ax.text(x+.12,y+h-.7,body,va='top',fontsize=9,linespacing=1.45)
def arrow(a,b):ax.add_patch(FancyArrowPatch(a,b,arrowstyle='-|>',mutation_scale=13,color='#596674',linewidth=1.2))
box(.3,6.6,9.4,1.3,'Same target, same endpoint, two assays','Source potency for every candidate + structure + provenance\nIC50 and Ki remain separate; identity / assay / document groups are retained','#1f4e70')
box(.3,3.7,4.3,2,'Fixed-reference transfer','Original 8 target references\nFixed analogue query graph\nPredict target potency; choose an analogue','#217c72')
box(5.2,3.7,4.5,2,'Reference measurement design','Fixed pool disjoint from fixed queries\nDiverse 4 + 4 selected measurements\nBudgetOracle: new = 1, repeat = 0','#956339')
arrow((2.5,6.6),(2.5,5.8));arrow((7.5,6.6),(7.5,5.8))
box(.3,1.5,9.4,1.3,'Freeze predictions and acquisition traces before query scoring','TaskView excludes query truth; a separate evaluator opens it afterwards.\nQuery / unknown / over-budget measurement requests are rejected and logged.','#1f4e70')
arrow((2.5,3.65),(2.5,2.85));arrow((7.5,3.65),(7.5,2.85))
ax.text(5,.7,'Selection loss + MAE + difference MSE',ha='center',weight='bold',fontsize=11)
ax.text(5,.23,'Weighted edges → equal tasks within group → equal groups',ha='center',fontsize=9)
arrow((5,1.45),(5,1.0));fig.tight_layout()
for ext in ['pdf','svg','png']:fig.savefig(O/('figure1_tracks.'+ext),dpi=220,bbox_inches='tight')
plt.close(fig)
ic=json.loads((ROOT/'reports/sar_transport_v1/ALL_METRICS.json').read_text())['models'];ki=json.loads((R/'ALL_METRICS.json').read_text());contr=json.loads((R/'PRIMARY_CONTRASTS.json').read_text())['comparisons']
fig=plt.figure(figsize=(7.1,6.4));gs=fig.add_gridspec(2,2,height_ratios=[1.6,1],hspace=.57,wspace=.8)
ax=fig.add_subplot(gs[0,0]);names=['source','affine','rbf','tanimoto','rf','fscap','actfound','no_condition','conditional','isotonic']
labels=['Source','Affine','Source RBF','Tanimoto','RF','FS-CAP*','ActFound*','Deep: no condition','Deep: conditional','Isotonic']
ax.scatter([ic[x]['L'] for x in names],range(10),c=['#217c72' if x=='tanimoto' else '#536575' for x in names],s=29)
ax.set_yticks(range(10),labels);ax.invert_yaxis();ax.set_xlabel('Selection loss (pIC50)');ax.set_title('A  IC50 fixed C8: 92 groups',loc='left',weight='bold',fontsize=10);ax.grid(axis='x',alpha=.2)
ax=fig.add_subplot(gs[0,1]);methods=['source','affine','rbf','tanimoto'];kk={x['model']:x for x in ki if x['track_id']=='fixed_c8_v1'}
ax.scatter([kk[x]['L'] for x in methods],range(4),c=['#536575','#536575','#536575','#217c72'],s=35)
ax.set_yticks(range(4),['Source','Affine','Source RBF','Tanimoto']);ax.invert_yaxis();ax.set_xlabel('Selection loss (pKi)');ax.set_title('B  Ki fixed C8: 50 groups',loc='left',weight='bold',fontsize=10);ax.grid(axis='x',alpha=.2)
ax=fig.add_subplot(gs[1,:]);vals=[contr[x]['difference_L'] for x in ['K1','K2']];lo=[contr[x]['family95_interval'][0] for x in ['K1','K2']];hi=[contr[x]['family95_interval'][1] for x in ['K1','K2']]
ax.errorbar(vals,[1,0],xerr=[np.array(vals)-lo,np.array(hi)-vals],fmt='o',color='#956339',capsize=4)
ax.axvline(0,color='#555555',linestyle='--',linewidth=1);ax.set_yticks([1,0],['K1: Tanimoto − source RBF\n(fixed C8)','K2: max-min − random mean\n(reference selection)']);ax.set_ylim(-.55,1.55)
ax.set_xlabel('Difference in selection loss (pKi); negative favors first method')
ax.set_title('C  Frozen Ki primary family: neither interval excludes zero',loc='left',fontsize=10,weight='bold');ax.grid(axis='x',alpha=.15)
fig.subplots_adjust(left=.28,right=.98,top=.96,bottom=.12)
for ext in ['pdf','svg','png']:fig.savefig(O/('figure2_results.'+ext),dpi=220,bbox_inches='tight')
plt.close(fig)
fig,ax=plt.subplots(figsize=(7,4))
styles={'maxmin':('#217c72','-'),'random_mean':('#536575','--'),'pairwise_kg':('#956339','-')}
for model in ['maxmin','facility','point_vopt','contrast_vopt','source_spread','pairwise_kg','random_mean']:
 rr=sorted([x for x in ki if x['track_id']=='acquisition_sentinel_v1' and x['model']==model],key=lambda x:x['budget'])
 color,ls=styles.get(model,(None,':'));ax.plot([x['budget'] for x in rr],[x['L'] for x in rr],label=model,marker='o',markersize=3,color=color,linestyle=ls)
ax.set(xlabel='Paid reference compounds (initial four included)',ylabel='Group-equal selection loss (pKi)',xticks=range(4,9));ax.legend(fontsize=8,ncol=2);ax.grid(alpha=.15);fig.tight_layout()
for ext in ['pdf','svg','png']:fig.savefig(O/('figureS1_budget_curves.'+ext),dpi=200,bbox_inches='tight')
plt.close(fig)
(R/'FIGURE_RENDER.json').write_text(json.dumps({'CPU_seconds':time.process_time()-start,'source':'Saved IC50 and frozen Ki metrics only','new_fits':0,'copied_publication_figures':False})+'\n')
