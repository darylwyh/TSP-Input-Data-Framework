'''
utility designed to log and report various statistics and performance metrics 
during the training process of a machine learning model. 
This function helps monitor the model's performance, gradient behaviors, 
and losses, essential for debugging and optimizing the training process

cost: A tensor containing the costs or losses for a batch of data.
grad_norms: A tuple containing the norms of gradients before and after clipping. This helps in monitoring gradient exploding or vanishing problems.
epoch: Current epoch number in training.
batch_id: Identifier for the current batch within the training dataset.
step: The overall step count in training, often equivalent to the number of batches processed.
log_likelihood: Tensor containing the log-likelihood values calculated during training, used to monitor the model’s performance in terms of probability estimation.
reinforce_loss: The loss calculated specifically from the REINFORCE algorithm part of the training if applicable.
bl_loss: The baseline loss, used when a baseline (like a critic in actor-critic methods) is employed to reduce variance in policy gradient estimates.
tb_logger: An instance of a TensorBoard logger used to log metrics to TensorBoard for visualization.
opts: A configuration object that might contain various settings and flags used throughout the training process.
'''
def log_values(cost, grad_norms, epoch, batch_id, step,
               log_likelihood, reinforce_loss, bl_loss, tb_logger, opts):
    # Computes the average cost from the cost tensor and converts it to a Python scalar with .item(). 
    # Provides a single representative value for the cost of the current batch, simplifying its logging and interpretation
    avg_cost = cost.mean().item()
    # Extracts the gradients norms into two separate variables for unclipped and clipped gradients. 
    # Monitoring both helps in understanding how gradient clipping affects the training dynamics
    grad_norms, grad_norms_clipped = grad_norms

    # Log values to screen
    print('epoch: {}, train_batch_id: {}, avg_cost: {}'.format(epoch, batch_id, avg_cost))
    # Logs the first component of the gradient norms (assuming multiple components might be present in some models, 
    # such as in networks with separate parts like actor and critic)
    print('grad_norm: {}, clipped: {}'.format(grad_norms[0], grad_norms_clipped[0]))

    # Log values to tensorboard
    if not opts.no_tensorboard:
        tb_logger.log_value('avg_cost', avg_cost, step)

        tb_logger.log_value('actor_loss', reinforce_loss.item(), step)
        tb_logger.log_value('nll', -log_likelihood.mean().item(), step)

        tb_logger.log_value('grad_norm', grad_norms[0], step)
        tb_logger.log_value('grad_norm_clipped', grad_norms_clipped[0], step)

        if opts.baseline == 'critic':
            tb_logger.log_value('critic_loss', bl_loss.item(), step)
            tb_logger.log_value('critic_grad_norm', grad_norms[1], step)
            tb_logger.log_value('critic_grad_norm_clipped', grad_norms_clipped[1], step)
