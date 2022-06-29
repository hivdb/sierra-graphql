/*

    Copyright (C) 2017 Stanford HIVDB team

    Sierra is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Sierra is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

package edu.stanford.hivdb.graphql;

import static edu.stanford.hivdb.graphql.ExtGraphQL.getPropertyViaMethod;
import static graphql.Scalars.GraphQLString;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

import org.apache.commons.lang3.tuple.Triple;

import edu.stanford.hivdb.drugs.DrugClass;
import edu.stanford.hivdb.mutations.MutationSet;
import edu.stanford.hivdb.mutations.MutationType;
import edu.stanford.hivdb.seqreads.SequenceReads;
import edu.stanford.hivdb.sequences.AlignedSequence;
import edu.stanford.hivdb.viruses.Gene;
import edu.stanford.hivdb.viruses.Virus;
import edu.stanford.hivdb.viruses.WithGene;
import graphql.execution.DataFetcherResult;
import graphql.schema.DataFetcher;
import graphql.schema.DataFetchingEnvironment;
import graphql.schema.GraphQLEnumType;
import graphql.schema.GraphQLFieldDefinition.Builder;
import graphql.schema.GraphQLList;
import graphql.schema.GraphQLTypeReference;

public class MutationSetDef {
	
	private enum mutsFilterOption {
		APOBEC, APOBEC_DRM,
		DRM, notDRM,
		SEQUENCED_ONLY,
		@Deprecated
		PI_DRM,
		@Deprecated
		NRTI_DRM,
		@Deprecated
		NNRTI_DRM,
		@Deprecated
		INSTI_DRM,
		SDRM, notSDRM,
		@Deprecated
		PI_SDRM,
		@Deprecated
		NRTI_SDRM,
		@Deprecated
		NNRTI_SDRM,
		@Deprecated
		INSTI_SDRM,
		DRP,
		TSM, notTSM,
		@Deprecated
		PI_TSM,
		@Deprecated
		NRTI_TSM,
		@Deprecated
		NNRTI_TSM,
		@Deprecated
		INSTI_TSM,
		@Deprecated
		GENE_PR,
		@Deprecated
		GENE_RT,
		@Deprecated
		GENE_IN,
		@Deprecated
		TYPE_MAJOR,
		@Deprecated
		TYPE_ACCESSORY,
		@Deprecated
		TYPE_NRTI,
		@Deprecated
		TYPE_NNRTI,
		@Deprecated
		TYPE_OTHER,
		INSERTION,
		DELETION,
		UNUSUAL,
		AMBIGUOUS,
		STOPCODON,
		CUSTOMLIST
	};

	private static GraphQLEnumType oMutationSetFilterOption =
		GraphQLEnumType.newEnum()
		.name("MutationSetFilterOption")
		.description("Filter option for mutation set.")
		.value(
			"APOBEC", mutsFilterOption.APOBEC,
			"List only mutations which are APOBEC-mediated " +
			"G-to-A hypermutation.")
		.value(
			"APOBEC_DRM", mutsFilterOption.APOBEC_DRM,
			"List only drug resistance mutations which are " +
			"APOBEC-mediated G-to-A hypermutation.")
		.value(
			"DRM", mutsFilterOption.DRM,
			"List only mutations which are drug resistance mutation (DRM).")
		.value(
			"DRP", mutsFilterOption.DRP,
			"List all mutations at DRM positions (no need to be DRMs).")
		.value(
			"notDRM", mutsFilterOption.notDRM,
			"List only mutations which are not drug resistance mutation (DRM).")
		.value(
			"SEQUENCED_ONLY", mutsFilterOption.SEQUENCED_ONLY,
			"Remove all unsequenced positions.")
		.value(
			"PI_DRM", mutsFilterOption.PI_DRM,
			"List only mutations which are PI DRM.",
			"Use combination of `drugClass=PI` and `filterOptions=DRM` instead.")
		.value(
			"NRTI_DRM", mutsFilterOption.NRTI_DRM,
			"List only mutations which are NRTI DRM.",
			"Use combination of `drugClass=NRTI` and `filterOptions=DRM` instead.")
		.value(
			"NNRTI_DRM", mutsFilterOption.NNRTI_DRM,
			"List only mutations which are NNRTI DRM.",
			"Use combination of `drugClass=NNRTI` and `filterOptions=DRM` instead.")
		.value(
			"INSTI_DRM", mutsFilterOption.INSTI_DRM,
			"List only mutations which are INSTI DRM.",
			"Use combination of `drugClass=INSTI` and `filterOptions=DRM` instead.")
		.value(
			"SDRM", mutsFilterOption.SDRM,
			"List only mutations which are surveillance drug resistance " +
			"mutations (SDRM).")
		.value(
			"notSDRM", mutsFilterOption.notSDRM,
			"List only mutations which are not sruveillance drug resistance " +
			"mutation (SDRM).")
		.value(
			"PI_SDRM", mutsFilterOption.PI_SDRM,
			"List only mutations which are PI SDRM.",
			"Use combination of `drugClass=PI` and `filterOptions=SDRM` instead.")
		.value(
			"NRTI_SDRM", mutsFilterOption.NRTI_SDRM,
			"List only mutations which are NRTI SDRM.",
			"Use combination of `drugClass=NRTI` and `filterOptions=SDRM` instead.")
		.value(
			"NNRTI_SDRM", mutsFilterOption.NNRTI_SDRM,
			"List only mutations which are NNRTI SDRM.",
			"Use combination of `drugClass=NNRTI` and `filterOptions=SDRM` instead.")
		.value(
			"INSTI_SDRM", mutsFilterOption.INSTI_SDRM,
			"List only mutations which are INSTI SDRM.",
			"Use combination of `drugClass=INSTI` and `filterOptions=SDRM` instead.")
		.value(
			"TSM", mutsFilterOption.TSM,
			"List only mutations which are treatment-selected mutations (TSM).")
		.value(
			"notTSM", mutsFilterOption.notTSM,
			"List only mutations which are not treatment-selected mutations (TSM).")
		.value(
			"PI_TSM", mutsFilterOption.PI_TSM,
			"List only mutations which are PI TSM.",
			"Use combination of `drugClass=PI` and `filterOptions=TSM` instead.")
		.value(
			"NRTI_TSM", mutsFilterOption.NRTI_TSM,
			"List only mutations which are NRTI TSM.",
			"Use combination of `drugClass=NRTI` and `filterOptions=TSM` instead.")
		.value(
			"NNRTI_TSM", mutsFilterOption.NNRTI_TSM,
			"List only mutations which are NNRTI TSM.",
			"Use combination of `drugClass=NNRTI` and `filterOptions=TSM` instead.")
		.value(
			"INSTI_TSM", mutsFilterOption.INSTI_TSM,
			"List only mutations which are INSTI TSM.",
			"Use combination of `drugClass=INSTI` and `filterOptions=TSM` instead.")
		.value("GENE_PR", mutsFilterOption.GENE_PR,
			"Use `includeGenes=PR instead.")
		.value("GENE_RT", mutsFilterOption.GENE_RT,
			"Use `includeGenes=RT instead.")
		.value("GENE_IN", mutsFilterOption.GENE_IN,
			"Use `includeGenes=IN instead.")
		.value("TYPE_MAJOR", mutsFilterOption.TYPE_MAJOR,
			"Use `mutationType=Major` instead.")
		.value("TYPE_ACCESSORY", mutsFilterOption.TYPE_ACCESSORY,
			"Use `mutationType=Accessory` instead.")
		.value("TYPE_NRTI", mutsFilterOption.TYPE_NRTI,
			"Use `mutationType=NRTI` instead.")
		.value("TYPE_NNRTI", mutsFilterOption.TYPE_NNRTI,
			"Use `mutationType=NNRTI` instead.")
		.value("TYPE_OTHER", mutsFilterOption.TYPE_OTHER,
			"Use `mutationType=Other` instead.")
		.value("INSERTION", mutsFilterOption.INSERTION)
		.value("DELETION", mutsFilterOption.DELETION)
		.value("UNUSUAL", mutsFilterOption.UNUSUAL)
		.value(
			"AMBIGUOUS", mutsFilterOption.AMBIGUOUS,
			"List all highly-ambiguous (HBDVN) mutations.")
		.value(
			"STOPCODON", mutsFilterOption.STOPCODON,
			"List only mutations with stop codon(s).")
		.value(
			"CUSTOMLIST", mutsFilterOption.CUSTOMLIST,
			"Accept a custom list of mutations and find the intersects.")
		.build();


	final public static <VirusT extends Virus<VirusT>> DataFetcher<DataFetcherResult<MutationSet<VirusT>>> newMutationSetDataFetcher(VirusT virusIns, String propertyName) {
		return env -> {
			Object src = env.getSource();
			@SuppressWarnings("unchecked")
			MutationSet<VirusT> mutations = (MutationSet<VirusT>) getPropertyViaMethod(propertyName, src, env);
			mutations = filterMutations(mutations, virusIns, env);

			return (
				new DataFetcherResult.Builder<MutationSet<VirusT>>()
				.data(mutations)
				.localContext(src)
				.build()
			);
		};
	}
	
	
	final private static <VirusT extends Virus<VirusT>> MutationSet<VirusT> filterMutations(MutationSet<VirusT> mutations, VirusT virusIns, DataFetchingEnvironment env) {
		List<?> filterOptions = env.getArgument("filterOptions");
		Collection<String> includeGenes = env.getArgument("includeGenes");
		String drugClassText = env.getArgument("drugClass");
		DrugClass<?> drugClass = virusIns.getDrugClass(drugClassText);
		MutationType<?> mutType = env.getArgument("mutationType");
		if (filterOptions == null) { filterOptions = new ArrayList<>(); }
		for (Object filterOption : filterOptions) {
			switch((mutsFilterOption) filterOption) {
			case APOBEC:
				mutations = mutations.getApobecMutations();
				break;
			case APOBEC_DRM:
				mutations = mutations.getApobecDRMs();
				break;
			case DRM:
				mutations = mutations.getDRMs();
				break;
			case DRP:
				mutations = mutations.getAtDRPMutations();
			case notDRM:
				mutations = mutations.subtractsBy(mutations.getDRMs());
				break;
			case SEQUENCED_ONLY:
				mutations = mutations.filterBy(mut -> !mut.isUnsequenced());
				break;
			case PI_DRM:
				mutations = mutations.getDRMs(virusIns.getDrugClass("PI"));
				break;
			case NRTI_DRM:
				mutations = mutations.getDRMs(virusIns.getDrugClass("NRTI"));
				break;
			case NNRTI_DRM:
				mutations = mutations.getDRMs(virusIns.getDrugClass("NNRTI"));
				break;
			case INSTI_DRM:
				mutations = mutations.getDRMs(virusIns.getDrugClass("INSTI"));
				break;
			case SDRM:
				mutations = mutations.getSDRMs();
				break;
			case notSDRM:
				mutations = mutations.subtractsBy(mutations.getSDRMs());
				break;
			case PI_SDRM:
				mutations = mutations.getSDRMs(virusIns.getDrugClass("PI"));
				break;
			case NRTI_SDRM:
				mutations = mutations.getSDRMs(virusIns.getDrugClass("NRTI"));
				break;
			case NNRTI_SDRM:
				mutations = mutations.getSDRMs(virusIns.getDrugClass("NNRTI"));
				break;
			case INSTI_SDRM:
				mutations = mutations.getSDRMs(virusIns.getDrugClass("INSTI"));
				break;
			case TSM:
				mutations = mutations.getTSMs();
				break;
			case notTSM:
				mutations = mutations.subtractsBy(mutations.getTSMs());
				break;
			case PI_TSM:
				mutations = mutations.getTSMs(virusIns.getDrugClass("PI"));
				break;
			case NRTI_TSM:
				mutations = mutations.getTSMs(virusIns.getDrugClass("NRTI"));
				break;
			case NNRTI_TSM:
				mutations = mutations.getTSMs(virusIns.getDrugClass("NNRTI"));
				break;
			case INSTI_TSM:
				mutations = mutations.getTSMs(virusIns.getDrugClass("INSTI"));
				break;
			case GENE_PR:
				// TODO: HIV2 support
				mutations = mutations.getGeneMutations(virusIns.getGene("HIV1PR"));
				break;
			case GENE_RT:
				mutations = mutations.getGeneMutations(virusIns.getGene("HIV1RT"));
				break;
			case GENE_IN:
				mutations = mutations.getGeneMutations(virusIns.getGene("HIV1IN"));
				break;
			case TYPE_MAJOR:
				mutations = mutations.getByMutType(virusIns.getMutationType("Major"));
				break;
			case TYPE_ACCESSORY:
				mutations = mutations.getByMutType(virusIns.getMutationType("Accessory"));
				break;
			case TYPE_NRTI:
				mutations = mutations.getByMutType(virusIns.getMutationType("NRTI"));
				break;
			case TYPE_NNRTI:
				mutations = mutations.getByMutType(virusIns.getMutationType("NNRTI"));
				break;
			case TYPE_OTHER:
				mutations = mutations.getByMutType(virusIns.getMutationType("Other"));
				break;
			case DELETION:
				mutations = mutations.getDeletions();
				break;
			case INSERTION:
				mutations = mutations.getInsertions();
				break;
			case UNUSUAL:
				mutations = mutations.getUnusualMutations();
				break;
			case AMBIGUOUS:
				mutations = mutations.getAmbiguousCodons();
				break;
			case STOPCODON:
				mutations = mutations.getStopCodons();
				break;
			case CUSTOMLIST:
				List<String> customList = env.getArgument("customList");
				Gene<VirusT> gene = null;
				if (WithGene.class.isInstance(env.getSource())) {
					WithGene<VirusT> source = env.getSource();
					gene = source.getGene();
				}
				MutationSet<VirusT> filterSet = virusIns.newMutationSet(gene, customList);
				mutations = mutations.intersectsWith(filterSet);
				break;
			default: break;
			}
		}
		if (includeGenes != null) {
			mutations = mutations.filterByNoSplit(
				mut -> includeGenes.contains(mut.getAbstractGene())
			);
		}
		if (drugClass != null) {
			mutations = mutations.filterBy(
				mut -> (
					mut.getDRMDrugClass() == drugClass ||
					mut.getSDRMDrugClass() == drugClass ||
					mut.getTSMDrugClass() == drugClass
				)
			);
		}
		if (mutType != null) {
			mutations = mutations.filterBy(mut -> mut.getPrimaryType() == mutType);
		}
		return mutations;
	}

	public static Builder newMutationSet(String virusName, Builder field, String name) {
		return field
			.name(name)
			.type(new GraphQLList(new GraphQLTypeReference("Mutation")))
			.argument(arg -> arg
				.name("filterOptions")
				.type(new GraphQLList(oMutationSetFilterOption))
				.description("List of filter options for the mutation set."))
			.argument(arg -> arg
				.name("includeGenes")
				.type(new GraphQLList(GeneDef.enumGene.get(virusName)))
				.defaultValue(Virus.getInstance(virusName).getAbstractGenes())
				.description("Specify a gene/genes for filtering the mutation set."))
			.argument(arg -> arg
				.name("drugClass")
				.type(DrugClassDef.enumDrugClass.get(virusName))
				.description("Specify a drug class for filtering the mutation set."))
			.argument(arg -> arg
				.name("mutationType")
				.type(MutationDef.enumMutationType.get(virusName))
				.description("Specify a mutation type for filtering the mutation set."))
			.argument(arg -> arg
				.name("customList")
				.type(new GraphQLList(GraphQLString))
				.description(
					"List of possible mutation strings that should be " +
					"included in this query if presented. Gene need to be " +
					"prepend if the gene is not able to be inferred from " +
					"the context."));
	}

	public static <VirusT extends Virus<VirusT>> MutationSet<VirusT> getMutationSetFromSource(Object src) {
		MutationSet<?> mutations;
		if (src instanceof AlignedSequence) {
			mutations = ((AlignedSequence<?>) src).getMutations(); 
		}
		else if (src instanceof SequenceReads) {
			mutations = ((SequenceReads<?>) src).getMutations();
		}
		else if (src instanceof Triple) {
			Object middle = ((Triple<?, ?, ?>) src).getMiddle();
			if (middle instanceof MutationSet) {
				mutations = (MutationSet<?>) middle;
			}
			else {
				throw new UnsupportedOperationException();
			}
		}
		else {
			throw new UnsupportedOperationException();
		}
		@SuppressWarnings("unchecked")
		MutationSet<VirusT> virusMutations = (MutationSet<VirusT>) mutations;
		return virusMutations;
	}
	
}
